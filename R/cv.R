
#' Cross-validating distributional anchor regression models
#'
#' @param x Object of class \code{"tranchor"}.
#' @param epochs Number of epochs.
#' @param folds Optional vector of length <number of samples> specifying the
#'     folds.
#' @param fold If \code{folds} is not specified, folds are computed from the
#'     anchors. See details.
#' @param ... Arguments passed to \code{\link[tranchor]{fit.tranchor}}.
#' @param xi Optional; can be used to overwrite \code{xi} used when building
#'     \code{x}.
#' @param verbose Whether to show a progress bar; defaults to \code{TRUE}.
#' @param print_training_info Whether to show the training steps in each fold;
#'     defaults to \code{FALSE}.
#'
#' @return List of negative log-likelihood contributions in- and out-of-sample,
#'     \code{folds} and \code{fold}.
#'
#' @exportS3Method cv tranchor
#'
cv.tranchor <- function(x, epochs, xi = NULL, folds = NULL, fold = 5,
                        verbose = TRUE, print_training_info = FALSE, ...) {
  call <- x$init_params$call
  data <- eval(call$data, envir = parent.frame())
  anchor <- eval(call$anchor, envir = parent.frame())
  if (!is.null(call$optimizer))
    call$optimizer <- eval(call$optimizer, envir = parent.frame())

  if (is.null(folds)) {
    folds <- .compute_folds(anchor, data, fold)
  }
  fold <- max(afolds <- sort(unique(as.numeric(as.factor(folds)))))

  out <- matrix(nrow = nrow(data), ncol = fold + 1)
  colnames(out) <- c(paste0("fold", afolds), "test")
  if (verbose && interactive())
    pb <- utils::txtProgressBar(max = fold, style = 3)
  for (tfold in afolds) {
    if (!is.null(x$init_params$call$optimizer))
      call$optimizer <- eval(x$init_params$call$optimizer,
                             envir = parent.frame())
    if (verbose && interactive())
      utils::setTxtProgressBar(pb, tfold)
    tridx <- which(folds != tfold)
    teidx <- which(folds == tfold)
    train <- data[tridx, ]
    test <- data[teidx, ]
    call$data <- train
    m <- eval(call)
    fit(m, epochs = epochs, verbose = print_training_info, ...)
    out[tridx, tfold] <- logLik(m, convert_fun = \(x) -x)
    out[teidx, fold + 1] <- logLik(m, newdata = test, convert_fun = \(x) -x)
  }

  structure(list(logLiki = out, folds = folds, fold = fold), class = "cvAN")
}

.compute_folds <- function(anchor, data, fold) {
  A <- .rm_int(model.matrix(anchor, data))
  if (ncol(A) == 1L) {
    if (length(unique(A)) > fold)
      folds <- cut(A, breaks = stats::quantile(A, probs = (0:fold)/fold),
                   include.lowest = TRUE)
    else folds <- factor(A)
  } else {
    folds <- as.factor(stats::cutree(stats::hclust(stats::dist(A)), k = fold))
  }
  as.numeric(folds)
}

.rm_int <- function(x) {
  if (is.matrix(x) && all(x[, 1] == 1))
    return(x[, -1, drop = FALSE])
  else
    x
}
