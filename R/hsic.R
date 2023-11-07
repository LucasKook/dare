
# tf_hsic <- function(x, y, kernel_x, kernel_y) {
#   m <- x$shape[[1]]
#   H <- tf$eye(m) - 1/m * tf$ones(c(m, m))
#   K <- kernel_x(x)
#   L <- kernel_y(y)
#   tf$linalg$diag_part(tf$matmul(L, tf$matmul(H, tf$matmul(K, H)))) / (m - 1)
# }

k_normal <- function(X, bw = 1) {
  sumX2 <- k_sum(X^2, axis = 2L, keepdims = TRUE)
  D2 <- sumX2 - 2.0 * tf$matmul(X, tf$transpose(X)) + tf$transpose(sumX2)
  D2 <- tf$clip_by_value(D2, 0, Inf)
  tf$exp(- D2 / (2 * bw^2))
}

indep_loss <- function(base_distribution, Amat, kx, ky, xi = 0) {

  if (is.character(base_distribution)) {
    bd <- deeptrafo:::get_bd(base_distribution)
  } else {
    bd <- base_distribution
  }

  Amat <- keras::k_constant(Amat)
  xi <- keras::k_constant(xi)
  kx <- k_constant(kx)
  ky <- k_constant(ky)

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

      m <- Amat$shape[[1]]
      H <- tf$eye(m) - 1/m * tf$ones(c(m, m))
      K <- k_normal(Amat, kx)
      # L <- k_normal(scores, ky)
      L <- k_normal(tfd_cdf(bd, trafo), ky)
      pen <- tf$linalg$diag_part(tf$matmul(L, tf$matmul(H, tf$matmul(K, H)))) / (m - 1)

      smpl <- k_flatten(tfd_cdf(bd, trafo))
      ecdf <- tfp$distributions$Empirical(smpl)
      unif <- (tfd_cdf(ecdf, smpl) - smpl)^2

      return(layer_add(list(10 * pen, unif))) # test xi = Inf
      # return(layer_add(list(neglogLik, xi * unif, xi * pen)))
    }
  )
}
