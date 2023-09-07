
devtools::load_all()
devtools::load_all("~/Dropbox/phd/projects/deeptrafo")

test_that("main works", {
  expect_no_error({
    set.seed(1)
    library("tram")

    dgp <- function(n = 100, g = \(x) qchisq(plogis(x), df = 3)) {
      A <- rnorm(n)
      X <- A + rnorm(n)
      Y <- g(X + A + rlogis(n))
      return(data.frame(Y = Y, X = X, A = A))
    }

    d <- dgp(n = tn <- 1e3)

    ### Sanity check
    tm <- Colr(Y ~ X, data = d, order = 10, support = range(d$Y))

    prm <- with(d, cbind(1, A) %*% solve(t(cbind(1, A)) %*% cbind(1, A)) %*% t(cbind(1, A)))
    m <- tranchor(Y ~ X, data = d, prm = prm, xi = 1e1, optimizer = optimizer_adam(0.1),
                  callbacks = list(callback_reduce_lr_on_plateau("loss")))
    fit(m, epochs = 1e3, batch_size = tn, validation_split = 0)

    logLik(tm)
    logLik.tranchor(m)

    unlist(coef(m))
    coef(tm)

    plot(residuals(m), residuals(tm), col = cut(d$A, breaks = quantile(d$A)))
  })
})

test_that("alias works", {
  expect_no_error({
    set.seed(1)
    library("tram")

    dgp <- function(n = 100, g = \(x) qchisq(plogis(x), df = 3)) {
      A <- rnorm(n)
      X <- A + rnorm(n)
      Y <- g(X + A + rlogis(n))
      return(data.frame(Y = Y, X = X, A = A))
    }

    d <- dgp(n = tn <- 1e3)

    ### Sanity check
    tm <- Colr(Y ~ X, data = d, order = 10, support = range(d$Y))

    prm <- with(d, cbind(1, A) %*% solve(t(cbind(1, A)) %*% cbind(1, A)) %*% t(cbind(1, A)))
    m <- ColrAN(Y ~ X, data = d, prm = prm, xi = 1e1, optimizer = optimizer_adam(0.1),
                callbacks = list(callback_reduce_lr_on_plateau("loss")))
    fit(m, epochs = 1e3, batch_size = tn, validation_split = 0)

    logLik(tm)
    logLik(m)

    unlist(coef(m))
    coef(tm)

    plot(residuals(m), residuals(tm), col = cut(d$A, breaks = quantile(d$A)))
  })
})
