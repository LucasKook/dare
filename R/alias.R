
ColrAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("ColrAN", class(ret))

  ret

}

CoxphAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("CoxphAN", class(ret))

  ret

}

LehmanAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gumbel", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("LehmanAN", class(ret))

  ret

}

BoxCoxAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("BoxCoxAN", class(ret))

  ret

}

PolrAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type == "ordered")

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("PolrAN", class(ret))

  ret

}

LmAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = 1L,
                                  response_type = response_type,
                                  y_basis_fun = eval_lin,
                                  y_basis_fun_lower = .empty_fun(eval_lin),
                                  y_basis_fun_prime = eval_lin_prime,
                                  basis = "shiftscale"),
    ...
) {

  stopifnot(response_type == "continuous")

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("LmAN", class(ret))

  ret

}

SurvregAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
    trafo_options = NULL,
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  if (response_type == "survival") {
    ybf <- function(y) eval_loglin(y[, 1])
    ybfl <- function(y) eval_loglin(y[, 1])
    ybfp <- function(y) eval_loglin_prime(y[, 1])
  } else {
    ybf <- eval_loglin
    ybfl <- .empty_fun(ybf)
    ybfp <- eval_loglin_prime
  }

  trafo_options <- trafo_control(
    order_bsp = 1L,
    response_type = response_type,
    y_basis_fun = ybf,
    y_basis_fun_lower = ybfl,
    y_basis_fun_prime = ybfp,
    basis = "shiftscale"
  )

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("SurvregAN", class(ret))

  ret

}

cotramAN <- function(
    formula, data,
    prm, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    ...
) {

  stopifnot(response_type == "count")

  tsupp <- range(data[[all.vars(formula)[1]]])
  trafo_options <- trafo_control(
    order_bsp = order,
    response_type = response_type,
    y_basis_fun = .get_eval_cotram(order, tsupp),
    y_basis_fun_lower = .get_eval_cotram_lower(order, tsupp),
    y_basis_fun_prime = .empty_fun(.get_eval_cotram(order, tsupp))
  )

  ret <- tranchor(formula = formula, data = data,
                  prm = prm, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  class(ret) <- c("cotramAN", class(ret))

  ret

}
