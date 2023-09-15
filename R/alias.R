
#' Continuous outcome logistic anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
ColrDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic",
    monitor_metrics = NULL, trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    ...
) {

  call <- match.call()
  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction, latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Cox proportional hazards anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
CoxphDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gompertz",
    monitor_metrics = NULL, trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    ...
) {

  call <- match.call()
  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction,
    latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Reverse-time proportional hazards anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
LehmannDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gumbel",
    monitor_metrics = NULL, trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    ...
) {

  call <- match.call()
  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction,
    latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Transformed-normal (Box-Cox-type) anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
BoxCoxDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "normal",
    monitor_metrics = NULL, trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    ...
) {

  call <- match.call()
  stopifnot(response_type %in% c("continuous", "survival"))

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction, latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Proportional odds logistic anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
PolrDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic",
    monitor_metrics = NULL, trafo_options = trafo_control(
      order_bsp = order, response_type = response_type),
    ...
) {

  call <- match.call()
  stopifnot(response_type == "ordered")

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction, latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Linear anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
LmDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "normal", monitor_metrics = NULL,
    ...
) {

  call <- match.call()
  stopifnot(response_type == "continuous")

  trop <- trafo_control(
    order_bsp = 1L,
    response_type = response_type,
    y_basis_fun = deeptrafo:::eval_lin,
    y_basis_fun_lower = deeptrafo:::.empty_fun(deeptrafo:::eval_lin),
    y_basis_fun_prime = deeptrafo:::eval_lin_prime,
    basis = "shiftscale")

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction,
    latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trop,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Parametric survival anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
SurvregDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "gompertz", monitor_metrics = NULL,
    trafo_options = NULL,
    ...
) {

  call <- match.call()
  stopifnot(response_type %in% c("continuous", "survival"))

  if (response_type == "survival") {
    ybf <- function(y) deeptrafo:::eval_loglin(y[, 1])
    ybfl <- function(y) deeptrafo:::eval_loglin(y[, 1])
    ybfp <- function(y) deeptrafo:::eval_loglin_prime(y[, 1])
  } else {
    ybf <- deeptrafo:::eval_loglin
    ybfl <- deeptrafo:::.empty_fun(ybf)
    ybfp <- deeptrafo:::eval_loglin_prime
  }

  trafo_options <- trafo_control(
    order_bsp = 1L,
    response_type = response_type,
    y_basis_fun = ybf,
    y_basis_fun_lower = ybfl,
    y_basis_fun_prime = ybfp,
    basis = "shiftscale"
  )

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction, latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}

#' Transformed-count anchor regression
#'
#' @inheritParams dare
#' @inheritDotParams dare
#' @return See \code{\link[dare]{dare}}
#'
#' @export
#'
cotramDA <- function(
    formula, data,
    anchor, xi = 0,
    response_type = get_response_type(data[[all.vars(formula)[1]]]),
    order = get_order(response_type, data[[all.vars(formula)[1]]]),
    addconst_interaction = 0, latent_distr = "logistic", monitor_metrics = NULL,
    ...
) {

  call <- match.call()
  stopifnot(response_type == "count")

  tsupp <- range(data[[all.vars(formula)[1]]])
  trafo_options <- trafo_control(
    order_bsp = order,
    response_type = response_type,
    y_basis_fun = deeptrafo:::.get_eval_cotram(order, tsupp),
    y_basis_fun_lower = deeptrafo:::.get_eval_cotram_lower(order, tsupp),
    y_basis_fun_prime = deeptrafo:::.empty_fun(
      deeptrafo:::.get_eval_cotram(order, tsupp)))

  ret <- dare(
    formula = formula, data = data,
    anchor = anchor, xi = xi,
    response_type = response_type, order = order,
    addconst_interaction = addconst_interaction, latent_distr = latent_distr,
    monitor_metrics = monitor_metrics, trafo_options = trafo_options,
    ... = ...)

  ret$init_params$call <- call
  ret

}
