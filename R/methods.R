
#' Log-likelihood for in- and out-of-sample model evaluation
#'
#' @exportS3Method logLik tranchor
#'
#' @param object Object of class \code{"tranchor"}
#' @param newdata Data.frame or list for evaluating on new data.
#' @param convert_fun See \code{\link[deeptrafo]{logLik.deeptrafo}}.
#' @param ... See \code{\link[deeptrafo]{logLik.deeptrafo}}.
#'
#' @export
logLik.tranchor <- function(
    object,
    newdata = NULL,
    convert_fun = function(x, ...) - sum(x, ...),
    ...
)
{

  if (object$init_params$is_atm && !is.null(newdata)) {
    lags <- deeptrafo:::fm_to_lag(object$init_params$lag_formula)
    newdata <- deeptrafo:::create_lags(rvar = object$init_params$response_varname,
                                       d_list = newdata,
                                       lags = lags)$data
  }

  if (is.null(newdata)) {
    y <- object$init_params$y
    y_pred <- deeptrafo:::fitted.deeptrafo(object, call_create_lags = FALSE, ... = ...)
  } else {
    y <- deeptrafo:::response(newdata[[object$init_params$response_varname]])
    y_pred <- deeptrafo:::fitted.deeptrafo(object, call_create_lags = FALSE,
                               newdata = newdata, ... = ...)
  }

  nll <- deeptrafo:::nll(object$init_params$latent_distr)
  convert_fun(nll(y, y_pred)$numpy())

}

.get_resid_fun <- function(base_distribution) {

  if (is.character(base_distribution)) {
    bd <- deeptrafo:::get_bd(base_distribution)
  } else {
    bd <- base_distribution
  }

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

      return(cleft * sc_left + exact * sc_exact + cright * sc_right + cint * sc_int)
    }
  )
}

#' Residuals for in- and out-of-sample model evaluation
#'
#' @param object Object of class \code{"tranchor"}
#' @param newdata Optional data.frame or list of new data to compute residuals on
#' @param convert_fun See \code{\link[deeptrafo]{residuals.deeptrafo}}
#' @param ... See \code{\link[deeptrafo]{residuals.deeptrafo}}
#'
#' @return Single-column matrix of score residuals
#'
#' @exportS3Method residuals tranchor
#' @export
residuals.tranchor <- function(
    object,
    newdata = NULL,
    convert_fun = function(x, ...) - identity(x),
    ...
)
{

  if (object$init_params$is_atm && !is.null(newdata)) {
    lags <- deeptrafo:::fm_to_lag(object$init_params$lag_formula)
    newdata <- deeptrafo:::create_lags(rvar = object$init_params$response_varname,
                                       d_list = newdata,
                                       lags = lags)$data
  }

  if (is.null(newdata)) {
    y <- object$init_params$y
    y_pred <- deeptrafo:::fitted.deeptrafo(object, call_create_lags = FALSE, ... = ...)
  } else {
    y <- deeptrafo:::response(newdata[[object$init_params$response_varname]])
    y_pred <- deeptrafo:::fitted.deeptrafo(object, call_create_lags = FALSE,
                               newdata = newdata, ... = ...)
  }

  res <- .get_resid_fun(object$init_params$latent_distr)
  convert_fun(res(y, y_pred)$numpy())

}

#' Fit distributional anchor regression models
#'
#' @param object Object of class \code{"tranchor"}.
#' @param epochs Number of epochs to train the model.
#' @param early_stopping_metric Metric monitored for early stopping. Defaults to
#'     \code{"loss"}.
#' @param callbacks List of optional callbacks.
#' @param ... See \code{\link[deepregression]{fit.deepregression}}
#'
#' @return See \code{\link[keras]{fit.keras.engine.training.Model}}
#'
#' @exportS3Method fit tranchor
#' @export
fit.tranchor <- function(
    object,
    epochs = 10,
    early_stopping_metric = "loss",
    callbacks = list(),
    ...
){
  deepregression::fit.deepregression(
    object, epochs = epochs, batch_size = nrow(object$init_params$y),
    shuffle = FALSE, early_stopping_metric = early_stopping_metric,
    callbacks = callbacks, validation_data = NULL, validation_split = 0,
    ...)
}
