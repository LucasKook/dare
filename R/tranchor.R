
#' Distributional anchor regression models
#'
#' @param formula Formula for specifying the regression of response on covariates
#' @param data Data.frame or list containing the variables occuring in \code{formula}
#'     and \code{anchor}
#' @param anchor Formula for specifying the anchor. The score residuals are
#'     projected onto the columns of the design matrix corresponding to
#'     \code{anchor} and \code{data}
#' @param xi Controls strength of causal regularization. \code{xi = 0} is the
#'     unpenalized model.
#' @param response_type See \code{\link[deeptrafo]{deeptrafo}}.
#' @param order See \code{\link[deeptrafo]{deeptrafo}}.
#' @param addconst_interaction See \code{\link[deeptrafo]{deeptrafo}}.
#' @param latent_distr See \code{\link[deeptrafo]{deeptrafo}}.
#' @param monitor_metrics See \code{\link[deeptrafo]{deeptrafo}}.
#' @param trafo_options See \code{\link[deeptrafo]{deeptrafo}}.
#' @param return_data See \code{\link[deeptrafo]{deeptrafo}}.
#' @param ... See \code{\link[deeptrafo]{deeptrafo}}.
#' @param loss Under development. Please use the default \code{loss = "anchor"}.
#' @param bw Under development.
#'
#' @return An untrained model of class \code{"dare"}.
#'
#' @importFrom mlt R
#' @importFrom Formula as.Formula
#' @importFrom stats model.matrix model.response model.frame dbeta as.formula
#'     fitted formula predict rmultinom logLik terms drop.terms
#' @importFrom keras layer_dense layer_add layer_concatenate
#' @export
dare <- function(
    formula,
    data,
    anchor,
    xi = 0,
    loss = c("anchor", "indep"),
    bw = c(1, 1),
    response_type = get_response_type(data[[all.vars(fml)[1]]]),
    order = get_order(response_type, data[[all.vars(fml)[1]]]),
    addconst_interaction = 0,
    latent_distr = "logistic",
    monitor_metrics = NULL,
    trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    return_data = FALSE,
    aggr = k_sum,
    ...
)
{

  call <- match.call()

  Amat <- stats::model.matrix(anchor, data)
  Q <- qr.Q(qr(Amat))
  prm <- tcrossprod(Q)

  # How many terms are in the formula
  fml <- Formula::as.Formula(formula)
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
        strsplit(deepregression::form2text(interacting), "+", fixed = TRUE)[[1]]
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
  y <- deeptrafo:::response(resp)

  # check for ATMs
  ftms <- attr(tms <- terms(list_of_formulas$h2), "term.labels")
  is_atm <- any(atps <- grepl("atplag", ftms))
  tlag_formula <- NULL
  if (is_atm) {
    # extract from lag formula the variables as simple sum and
    # layers for additional transformation
    tlag_formula <- paste0(grep("atplag", ftms, value = TRUE), collapse = "+")
    lags <- deeptrafo:::create_lags(rvar = rvar, d_list = data, atplags = tlag_formula)
    data <- lags$data

    resp <- data[[rvar]] # creating lags reduces data set size
    y <- deeptrafo:::response(resp)

    tlag_formula <- lags$fm
    list_of_formulas$yterms <- as.formula(
      paste0(deepregression::form2text(list_of_formulas$yterms), " + ", tlag_formula))
    if (length(ftms) > length(which(atps)))
      list_of_formulas$h2 <- stats::drop.terms(tms, which(atps))
    else
      list_of_formulas$h2 <- ~1
  }

  # define how to get a trafo model from predictor
  from_pred_to_trafo_fun <- deeptrafo:::from_preds_to_trafo(
    atm_toplayer = trafo_options$atm_toplayer, const_ia = addconst_interaction)

  atm_lag_processor <- deeptrafo:::atm_lag_processor_factory(rvar)

  trafo_processor <- list(
    bsfun = deeptrafo:::basis_processor, bsfunl = deeptrafo:::basis_processor_lower,
    bspfun = deeptrafo:::basisprime_processor, ia = deeptrafo:::ia_processor,
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
  loss <- match.arg(loss)
  if (loss == "anchor")
    tloss <- dare_loss(latent_distr, prm, xi)
  else if (loss == "indep")
    tloss <- indep_loss(latent_distr, Amat, bw[1], bw[2], xi, aggr = aggr)

  snwb <- list(subnetwork_init)[rep(1, length(list_of_formulas))]
  snwb[[which(names(list_of_formulas) == "h1pred")]] <-
    deeptrafo:::h1_init(yterms = which(names(list_of_formulas) == "yterms"),
            h1pred = which(names(list_of_formulas) == "h1pred"),
            add_const_positiv = addconst_interaction)
  snwb[[which(names(list_of_formulas) == "yterms")]] <- function(...)
    return(NULL)

  args <- c(list(
    y = y, family = latent_distr, data = data, list_of_formulas = list_of_formulas,
    subnetwork_builder = snwb, from_preds_to_output = from_pred_to_trafo_fun,
    loss = tloss, monitor_metrics = monitor_metrics,
    additional_processor = additional_processor, engine = "tf"), dots)

  ret <- suppressWarnings(do.call("deepregression", args))

  ret$init_params$is_atm <- is_atm
  ret$init_params$lag_formula <- tlag_formula
  ret$init_params$formula <- formula
  ret$init_params$trafo_options <- trafo_options
  ret$init_params$response_varname <- rvar
  ret$init_params$response_type <- response_type
  ret$init_params$response <- resp
  ret$init_params$prepare_y_valdata <- deeptrafo:::response
  ret$init_params$data <- if (return_data) data else NULL
  ret$init_params$latent_distr <- latent_distr
  ret$init_params$call <- call

  class(ret) <- c("dare", "deeptrafo", "deepregression")
  ret

}

#' @import deepregression
#' @import deeptrafo
#' @import tfprobability
#' @import keras
#' @import tensorflow
dare_loss <- function(base_distribution, prm, xi = 0) {

  if (is.character(base_distribution)) {
    bd <- deeptrafo:::get_bd(base_distribution)
  } else {
    bd <- base_distribution
  }

  prm <- keras::k_constant(prm)
  xi <- keras::k_constant(xi)

  return(
    function(y_true, y_pred) {

      cleft <- deepregression::tf_stride_cols(y_true, 1L)
      exact <- deepregression::tf_stride_cols(y_true, 2L)
      cright <- deepregression::tf_stride_cols(y_true, 3L)
      cint <- deepregression::tf_stride_cols(y_true, 4L)

      trafo <- keras::layer_add(list(deepregression::tf_stride_cols(y_pred, 1L), # Shift in 1
                              deepregression::tf_stride_cols(y_pred, 2L))) # Upper in 2
      trafo_lwr <- keras::layer_add(list(deepregression::tf_stride_cols(y_pred, 1L),
                                  deepregression::tf_stride_cols(y_pred, 3L))) # Lower in 3
      trafo_prime <- tf$math$log(tf$clip_by_value(deepregression::tf_stride_cols(y_pred, 4L),
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

      sc_exact <- dd / tf$clip_by_value(tfd_prob(bd, trafo), 1e-6, Inf)
      sc_left <- tfd_prob(bd, trafo) / tf$clip_by_value(tfd_cdf(bd, trafo), 1e-6, 1)
      sc_right <- - tfd_prob(bd, trafo_lwr) / tf$clip_by_value(1 - tfd_cdf(bd, trafo_lwr), 1e-6, 1)
      sc_int <- (tfd_prob(bd, trafo) - tfd_prob(bd, trafo_lwr)) /
        tf$clip_by_value(tfd_cdf(bd, trafo) - tfd_cdf(bd, trafo_lwr), 1e-16, 1)

      scores <- (cleft * sc_left + exact * sc_exact + cright * sc_right + cint * sc_int)
      pen <- tf$linalg$matmul(prm, scores)^2

      # return(pen) # test xi = Inf
      return(layer_add(list(neglogLik, xi * pen)))
    }
  )
}

