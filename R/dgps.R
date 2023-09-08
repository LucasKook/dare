
#' @export
simple_dgp <- function(n = 100, g = \(x) qchisq(plogis(x), df = 3)) {
  A <- rnorm(n)
  X <- A + rnorm(n)
  Y <- g(X + A + rlogis(n))
  return(data.frame(Y = Y, X = X, A = A))
}
