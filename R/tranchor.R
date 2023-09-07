
tranchor <- function(
    formula,
    data,
    prm,
    xi = 0,
    response_type = get_response_type(data[[all.vars(fml)[1]]]),
    order = get_order(response_type, data[[all.vars(fml)[1]]]),
    addconst_interaction = 0,
    latent_distr = "logistic",
    monitor_metrics = NULL,
    trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    return_data = FALSE,
    ...
)
{

  call <- match.call()

  # How many terms are in the formula
  fml <- as.Formula(formula)
  ninteracting <- length(attr(fml, "lhs"))
  nterms <- length(attr(fml, "rhs"))

  # Name of the response variable
  rvar <- all.vars(formula)[1]

  # Placeholder Intercept
  int <- 1

  # Set up formulas for basis
  if (ninteracting > 1L) {
    interacting <- formula(fml, lhs = 2L, rhs = 0L)[[2]]
    h1_form <- paste0(
      "~ 0 + ", paste(paste0("ia(", c(int, trimws(
        strsplit(form2text(interacting), "+", fixed = TRUE)[[1]]
      )), ")"), collapse=" + ")
    )
  } else {
    h1_form <- paste0(
      "~ -1 + ", paste(paste0("ia(", int, ")"), collapse=" + ")
    )
  }

  # List of formulas
  list_of_formulas <- list(
    yterms = as.formula(paste0("~ -1 + bsfun(", rvar, ") + bsfunl(", rvar,
                               ") + bspfun(", rvar, ")")),
    h1pred = as.formula(h1_form),
    h2 = if (nterms >= 1L) formula(fml, lhs = 0, rhs = 1L) else NULL,
    shared = if (nterms == 2L) formula(fml, lhs = 0, rhs = 2L) else NULL
  )

  # Remove NULL formulae
  list_of_formulas[sapply(list_of_formulas, is.null)] <- NULL

  # Extract response variable
  resp <- data[[rvar]]
  y <- response(resp)

  # check for ATMs
  ftms <- attr(tms <- terms(list_of_formulas$h2), "term.labels")
  is_atm <- any(atps <- grepl("atplag", ftms))
  tlag_formula <- NULL
  if (is_atm) {
    # extract from lag formula the variables as simple sum and
    # layers for additional transformation
    tlag_formula <- paste0(grep("atplag", ftms, value = TRUE), collapse = "+")
    lags <- create_lags(rvar = rvar, d_list = data, atplags = tlag_formula)
    data <- lags$data

    resp <- data[[rvar]] # creating lags reduces data set size
    y <- response(resp)

    tlag_formula <- lags$fm
    list_of_formulas$yterms <- as.formula(
      paste0(form2text(list_of_formulas$yterms), " + ", tlag_formula))
    if (length(ftms) > length(which(atps)))
      list_of_formulas$h2 <- drop.terms(tms, which(atps))
    else
      list_of_formulas$h2 <- ~1
  }

  # define how to get a trafo model from predictor
  from_pred_to_trafo_fun <- from_preds_to_trafo(
    atm_toplayer = trafo_options$atm_toplayer, const_ia = addconst_interaction)

  atm_lag_processor <- atm_lag_processor_factory(rvar)

  trafo_processor <- list(
    bsfun = basis_processor, bsfunl = basis_processor_lower,
    bspfun = basisprime_processor, ia = ia_processor,
    atplag = atm_lag_processor)

  dots <- list(...)

  if (is.null(dots$additional_processor)) {

    additional_processor <- trafo_processor

  } else{

    additional_processor <- c(list(...)$additional_processor, trafo_processor)
    dots$additional_processor <- NULL

  }

  attr(additional_processor, "controls") <- trafo_options

  # Loss function
  # tloss <- get_loss(response_type, latent_distr)
  tloss <- tranchor_loss(latent_distr, prm, xi)

  snwb <- list(subnetwork_init)[rep(1, length(list_of_formulas))]
  snwb[[which(names(list_of_formulas) == "h1pred")]] <-
    h1_init(yterms = which(names(list_of_formulas) == "yterms"),
            h1pred = which(names(list_of_formulas) == "h1pred"),
            add_const_positiv = addconst_interaction)
  snwb[[which(names(list_of_formulas) == "yterms")]] <- function(...)
    return(NULL)

  args <- c(list(
    y = y, family = latent_distr, data = data, list_of_formulas = list_of_formulas,
    subnetwork_builder = snwb, from_preds_to_output = from_pred_to_trafo_fun,
    loss = tloss, monitor_metrics = monitor_metrics,
    additional_processor = additional_processor), dots)

  ret <- suppressWarnings(do.call("deepregression", args))

  ret$init_params$is_atm <- is_atm
  ret$init_params$lag_formula <- tlag_formula
  ret$init_params$formula <- formula
  ret$init_params$trafo_options <- trafo_options
  ret$init_params$response_varname <- rvar
  ret$init_params$response_type <- response_type
  ret$init_params$response <- resp
  ret$init_params$prepare_y_valdata <- response
  ret$init_params$data <- if (return_data) data else NULL
  ret$init_params$call <- call

  class(ret) <- c("deeptrafo", "deepregression")
  ret

}
tranchor_loss <- function(base_distribution, prm, xi = 0) {

  if (is.character(base_distribution)) {
    bd <- deeptrafo:::get_bd(base_distribution)
  } else {
    bd <- base_distribution
  }

  prm <- k_constant(prm)
  xi <- k_constant(xi)

  return(
    function(y_true, y_pred) {

      cleft <- tf_stride_cols(y_true, 1L)
      exact <- tf_stride_cols(y_true, 2L)
      cright <- tf_stride_cols(y_true, 3L)
      cint <- tf_stride_cols(y_true, 4L)

      trafo <- layer_add(list(tf_stride_cols(y_pred, 1L), # Shift in 1
                              tf_stride_cols(y_pred, 2L))) # Upper in 2
      trafo_lwr <- layer_add(list(tf_stride_cols(y_pred, 1L),
                                  tf_stride_cols(y_pred, 3L))) # Lower in 3
      trafo_prime <- tf$math$log(tf$clip_by_value(tf_stride_cols(y_pred, 4L),
                                                  1e-8, Inf)) # Prime in 4

      ll_exact <- tfd_log_prob(bd, trafo) + trafo_prime
      ll_left <- tf$math$log(tf$clip_by_value(tfd_cdf(bd, trafo), 1e-16, 1))
      ll_right <- tf$math$log(tf$clip_by_value(1 - tfd_cdf(bd, trafo_lwr), 1e-16, 1))
      ll_int <- tf$math$log(tf$clip_by_value(tfd_cdf(bd, trafo) - tfd_cdf(bd, trafo_lwr), 1e-16, 1))

      neglogLik <- - (cleft * ll_left + exact * ll_exact + cright * ll_right + cint * ll_int)

      tape <- \() NULL
      with(tf$GradientTape() %as% tape, {
        tape$watch(trafo)
        dd2d <- tfd_prob(bd, trafo)
      })
      dd <- tape$gradient(dd2d, trafo)

      sc_exact <- dd / tf$clip_by_value(tfd_prob(bd, trafo), 1e-6, 20)
      sc_left <- tfd_prob(bd, trafo) / tf$clip_by_value(tfd_cdf(bd, trafo), 1e-6, 1)
      sc_right <- - tfd_prob(bd, trafo_lwr) / tf$clip_by_value(1 - tfd_cdf(bd, trafo_lwr), 1e-6, 1)
      sc_int <- (tfd_prob(bd, trafo) - tfd_prob(bd, trafo_lwr)) /
        tf$clip_by_value(tfd_cdf(bd, trafo) - tfd_cdf(bd, trafo_lwr), 1e-16, 1)

      scores <- (cleft * sc_left + exact * sc_exact + cright * sc_right + cint * sc_int)
      pen <- xi * tf$linalg$matmul(prm, scores)^2

      return(layer_add(list(neglogLik, pen)))
    }
  )
}

