
#' Continuous outcome logistic anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
ColrAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Cox proportional hazards anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
CoxphAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Reverse-time proportional hazards anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
LehmannAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gumbel", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Transformed-normal (Box-Cox-type) anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
BoxCoxAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- tranchor(formula = formula, data = data,
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Proportional odds logistic anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
PolrAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    trafo_options = trafo_control(order_bsp = order, response_type = response_type),
    ...
) {

  stopifnot(response_type == "ordered")

  ret <- tranchor(formula = formula, data = data,
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Linear anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
LmAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
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
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Parametric survival anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
SurvregAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
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
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}

#' Transformed-count anchor regression
#'
#' @inheritParams tranchor
#' @inheritDotParams tranchor
#' @return See \code{\link[tranchor]{tranchor}}
#'
#' @export
#'
cotramAN <- function(
    formula, data,
    anchor, xi = 0,
    response_type = tranchor:::get_response_type(data[[all.vars(formula)[1]]]),
    order = tranchor:::get_order(response_type, data[[all.vars(formula)[1]]]),
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
                  anchor = anchor, xi = xi,
                  response_type = response_type, order = order,
                  addconst_interaction = addconst_interaction, latent_distr = latent_distr,
                  monitor_metrics = monitor_metrics, trafo_options = trafo_options,
                  ... = ...)

  ret

}
