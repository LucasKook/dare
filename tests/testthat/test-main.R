
test_that("main works", {
  expect_no_error({
    set.seed(1)
    library("tram")

    d <- simple_dgp(n = tn <- 1e3)

    ### Sanity check
    tm <- Colr(Y ~ X, data = d, order = 10, support = range(d$Y))

    m <- tranchor(Y ~ X, data = d, anchor = ~ A, xi = 1e4,
                  optimizer = optimizer_adam(0.05),
                  callbacks = list(callback_reduce_lr_on_plateau("loss")))
    fit.tranchor(m, epochs = 1e3)

    logLik(tm)
    logLik.tranchor(m)
    logLik.tranchor(m, newdata = d[1:5, ])

    cbind(deeptrafo = c(unlist(coef(m, which = "inter")) + unlist(coef(m))[2],
                        unlist(coef(m)[-2])),
          tram = coef(tm, with_baseline = TRUE))

    # plot(residuals.tranchor(m), residuals(tm), col = cut(d$A, breaks = quantile(d$A)))
    # plot(residuals.tranchor(m) ~ A, data = d)
    # abline(lm(residuals.tranchor(m) ~ A, data = d))
    # plot(residuals(tm) ~ d$A)
    # abline(lm(residuals(tm) ~ A, data = d))
  })
})

test_that("alias works", {
  expect_no_error({
    set.seed(1)
    library("tram")

    d <- simple_dgp(n = tn <- 1e3)

    ### Sanity check
    tm <- Colr(Y ~ X, data = d, order = 10, support = range(d$Y))

    m <- ColrAN(Y ~ X, data = d, anchor = ~ A, xi = 1e1,
                optimizer = optimizer_adam(0.1),
                callbacks = list(callback_reduce_lr_on_plateau("loss")))
    fit.tranchor(m, epochs = 1e3)

    logLik(tm)
    logLik(m)

    unlist(coef(m))
    coef(tm)

    # plot(residuals.tranchor(m), residuals(tm), col = cut(d$A, breaks = quantile(d$A)))
    # plot(residuals.tranchor(m) ~ A, data = d)
    # abline(lm(residuals.tranchor(m) ~ A, data = d))
  })
})

test_that("cv works with all kinds of anchors", {
  expect_no_error({
    set.seed(1)
    d <- simple_dgp(n = tn <- 1e3)
    m <- ColrAN(Y ~ X, data = d, anchor = ~ A, xi = 1)
    out <- cv(m, epochs = 1, xi = 0, fold = 2)
    d$dA <- factor(sample(-3:-1, tn, TRUE))
    md <- ColrAN(Y ~ X, data = d, anchor = ~ dA, xi = 1)
    outd <- cv(md, epochs = 1, xi = 0, fold = 2)
    mv <- ColrAN(Y ~ X, data = d, anchor = ~ A + dA, xi = 1)
    outv <- cv(mv, epochs = 1, xi = 0, fold = 2)
  })
})
