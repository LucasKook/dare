# README demo
# LK 2023

devtools::load_all()

# Training data -----------------------------------------------------------

set.seed(42)
n <- 1e3
A <- rnorm(n)
H <- A + rnorm(n)
X <- H + A + rnorm(n)
Y <- qchisq(plogis(H + X + A + rlogis(n)), df = 3)
train <- data.frame(Y = Y, X = X, A = A)

# Penalized model ---------------------------------------------------------

m <- ColrAN(Y ~ X, data = train, anchor = ~ A, xi = 10,
            optimizer = optimizer_adam(0.1))
fit(m, epochs = 1e4)
unlist(coef(m))
###           X (Intercept)
###   -1.309156   -1.294785
logLik(m)
### [1] -2175.023

# Unpenalized model -------------------------------------------------------

m0 <- ColrAN(Y ~ X, data = train, anchor = ~ A, xi = 0,
             optimizer = optimizer_adam(0.1))
fit(m0, epochs = 1e4)
unlist(coef(m0))
###           X (Intercept)
###   -1.216760   -1.319929
logLik(m0)
### [1] -2162.867

# Test data ---------------------------------------------------------------

A <- 1 + rnorm(n)
H <- A + rnorm(n)
X <- H + A + rnorm(n)
Y <- qchisq(plogis(H + X + A + rlogis(n)), df = 3)
test <- data.frame(Y = Y, X = X, A = A)

# Out-of-sample -----------------------------------------------------------

logLik(m0, newdata = test)
### [1] -2692.289
logLik(m, newdata = test)
### [1] -2659.036

# Plot --------------------------------------------------------------------

plot(ecdf(logLik(m0, newdata = test, convert_fun = identity)), xlim = c(0, 6),
     cex = 0.1, xlab = "negative log-likelihood contribution")
plot(ecdf(logLik(m, newdata = test, convert_fun = identity)), add = TRUE,
     xlim = c(0, 6), cex = 0.1, col = 2)
legend("bottomright", c("unpenalized", "penalized"), col = c(1, 2), lwd = 2)

cutA <- ordered(cut(test$A, breaks = quantile(test$A), include.lowest = TRUE))
plot(residuals(m, newdata = test), residuals(m0, newdata = test),
     col = cutA)
legend("bottomright", levels(cutA), lwd = 2, col = 1:length(levels(cutA)),
       title = "A in test data")

# Cross-validation --------------------------------------------------------

cvd <- cv(m, epochs = 1e4, xi = 100)
sum(cvd$logLiki[, "test"])
boxplot(cvd$logLiki)
boxplot(cvd$logLiki[, "test"] ~ cvd$folds)
