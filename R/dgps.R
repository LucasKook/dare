
#' Simple data-generating process
#'
#' @param n sample size
#' @param g inverse transformation function
#' @param shift shift in the anchor
#'
#' @return data.frame with Y, X, A
#'
#' @importFrom stats plogis qchisq rlogis rnorm
#' @export
#'
simple_dgp <- function(n = 100, g = \(x) qchisq(plogis(x), df = 3), shift = 0) {
  A <- shift + rnorm(n)
  H <- A + rnorm(n)
  X <- H + A + rnorm(n)
  Y <- g(H + X + A + rlogis(n))
  return(data.frame(Y = Y, X = X, A = A))
}
