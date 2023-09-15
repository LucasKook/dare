
#' Simple data-generating process
#'
#' @param n sample size
#' @param g inverse transformation function
#' @param shift shift in the anchor
#' @param iv whether the graphical iv assumptions should be met
#'
#' @return data.frame with Y, X, A, H
#'
#' @importFrom stats plogis qchisq rlogis rnorm
#' @export
#'
simple_dgp <- function(n = 100, g = \(x) qchisq(plogis(x), df = 3),
                       shift = 0, iv = FALSE) {
  A <- shift + rnorm(n)
  H <- (1 - iv) * A + rnorm(n)
  X <- H + A + rnorm(n)
  Y <- g(H + X + (1 - iv) * A + rlogis(n))
  data.frame(Y = Y, X = X, A = A, H = H)
}
