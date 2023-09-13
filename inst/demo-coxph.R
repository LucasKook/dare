# Coxph anchor example
# LK 2023

# Dependencies ------------------------------------------------------------

library("tranchor")
library("tram")
library("survival")

# Data --------------------------------------------------------------------

data("GBSG2", package = "TH.data")
GBSG2$surv <- with(GBSG2, Surv(time, cens))

train <- GBSG2[GBSG2$age < quantile(GBSG2$age, probs = 0.5), ]
test <- GBSG2[GBSG2$age >= quantile(GBSG2$age, probs = 0.5), ]

# Unpenalized Coxph -------------------------------------------------------

m <- BoxCoxAN(surv ~ horTh + pnodes, data = train, anchor = ~ age, xi = 0,
              optimizer = optimizer_adam(0.1))
fit(m, epochs = 1e4)

# Anchor ------------------------------------------------------------------

ma <- BoxCoxAN(surv ~ horTh + pnodes, data = train, anchor = ~ age, xi = 1.5,
               optimizer = optimizer_adam(0.1))
fit(ma, epochs = 1e4)

# Visualize ---------------------------------------------------------------

opar <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
### No penalization -> residuals are correlated with anchor age
plot(residuals(m) ~ age, data = train, xlim = range(GBSG2$age), main = "unpenalized")
abline(lm(residuals(m) ~ age, data = train))

points(residuals(m, newdata = test) ~ age, data = test, col = 2)
abline(lm(residuals(m, newdata = test) ~ age, data = test), col = 2)

### With penalization -> residuals are more uncorrelated with anchor age
plot(residuals(ma) ~ age, data = train, xlim = range(GBSG2$age), main = "penalized")
abline(lm(residuals(ma) ~ age, data = train))

points(residuals(ma, newdata = test) ~ age, data = test, col = 2)
abline(lm(residuals(ma, newdata = test) ~ age, data = test), col = 2)
par(opar)

# GOF ---------------------------------------------------------------------

logLik(m)
logLik(ma)

# Out-of-sample prediction ------------------------------------------------

logLik(m, newdata = test)
logLik(ma, newdata = test)

# Individual NLL contributions --------------------------------------------

plot(ecdf(logLik(m, newdata = test, convert_fun = identity)))
lines(ecdf(logLik(ma, newdata = test, convert_fun = identity)), col = 2)
legend("topleft", c("unpenalized", "penalized"), col = 1:2, lwd = 2)

# Coefficients ------------------------------------------------------------

unlist(coef(m))
unlist(coef(ma))
