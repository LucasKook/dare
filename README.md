<!-- badges: start -->
  [![R-CMD-check](https://github.com/LucasKook/dare/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/dare/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Distributional anchor regression in R

The `dare` package implements distributional anchor regression [1] using a
(stochastic) gradient descent optimizer. Anchor regression [2] is a method to
robustify predictions under distribution shift induced by so-called anchors.
A small illustration is given below.

# Using package `dare`

Consider the following data generated from a structural causal model in which
`A` is the anchor, `X` a predictor, `H` a hidden confounder, and `Y` the 
response:
```r
set.seed(42)
n <- 1e3
A <- rnorm(n)
H <- A + rnorm(n)
X <- H + A + rnorm(n)
Y <- qchisq(plogis(H + X + A + rlogis(n)), df = 3)
train <- data.frame(Y = Y, X = X, A = A)
```

Conditional on `H + X + A = 0`, this generates responses with a Chi-square
distribution with three degrees of freedom. We set up and fit continuous outcome
logistic anchor regression below with regularization parameter `xi = 10`:
```r
m <- ColrDA(Y ~ X, data = train, anchor = ~ A, xi = 10, 
            optimizer = optimizer_adam(0.1))
fit(m, epochs = 1e4)
unlist(coef(m))
###           X (Intercept) 
###   -1.309156   -1.294785 
logLik(m)
### [1] -2175.023
```

The unpenalized model (`xi = 0`) yields different coefficient estimates and
a better in-sample log-likelihood:
```r
m0 <- ColrDA(Y ~ X, data = train, anchor = ~ A, xi = 0, 
             optimizer = optimizer_adam(0.1))
fit(m, epochs = 1e4)
unlist(coef(m0))
###           X (Intercept) 
###   -1.216760   -1.319929 
logLik(m0)
### [1] -2162.867
```

We now generate some test data from a different environment, in which the
distribution of `A` is shifted, which in turn induces distribution shift in the
other covariates:
```r
A <- 1 + rnorm(n)
H <- A + rnorm(n)
X <- H + A + rnorm(n)
Y <- qchisq(plogis(H + X + A + rlogis(n)), df = 3)
test <- data.frame(Y = Y, X = X, A = A)
```

Below we show the performance on the shifted out-of-sample data. De-correlating
the score residuals from the anchors (induced by causal regularization) leads
to a better out-of-sample log-likelihood on shifted data;
```r
logLik(m0, newdata = test)
### [1] -2692.289
logLik(m, newdata = test)
### [1] -2659.036
```

Cross-validation provides one way to choose the strength of regularization
`xi`. In `dare`, cross-validation is implemented as follows: If there is
a single anchor variable and the variable is a factor, "leave-one-level-out"
cross validation is performed. If the single anchor variable is numeric, it
is cut into intervals with boundaries given by quantiles. For multivariable
anchors, the vector of `folds` has to be provided per observation. The code
below shows how cross-validation can be used in our example above:
```r
cvd <- cv(m, epochs = 1e4, xi = 100)
-sum(cvd$logLiki[, "test"])
### [1] -2865.85
```
The output of `cv()` contains a matrix `logLiki` of log-likelihood contributions 
for each fold and the test data.

The code for reproducing this example (together with two plots visualizing
out-of-sample log-likelihood contributions and score residuals for the two
models) can be found in `./inst/demo-readme.R`. Another demo for Cox 
proportional hazard anchor regression can be found in `./inst/demo-coxph.R`.

# Implemented model classes

For general outcome types and error distributions, the `dare()` function
can be used. However, for archetypal model classes we provide the following
alias (using the name of function corresponding to a model class and appending
`DA` for "Distributional Anchor"):

| **Function alias**  | **Corresponding model**    |
|---------------------|----------------------------|
| `BoxCoxDA()`        | `tram::BoxCox()`           | 
| `ColrDA()`          | `tram::Colr()`             |
| `cotramDA()`        | `cotram::cotram()`         |
| `CoxphDA()`         | `tram::Coxph()`            |
| `LehmannDA()`       | `tram::Lehmann()`          |
| `LmDA()`            | `tram::Lm()`               |
| `PolrDA()`          | `tram::Polr()`             |
| `SurvregDA()`       | `tram::Survreg()`          |

# Installing package `dare`

The package can be installed via:
```r
remotes::install_github("LucasKook/dare")
```

The implementation is build upon package `deeptrafo` [3] which is available on
[CRAN](https://CRAN.R-project.org/package=deeptrafo). For trouble-shooting
`Python`, `tensorflow`, and `tfprobability` installations, please see
[the `deepregression` documentation](https://github.com/neural-structured-additive-learning/deepregression#troubleshooting).

# References

[1] Kook, L., Sick, B., & Bühlmann, P. (2022). Distributional anchor regression. Statistics and Computing, 32(3), 39. [doi:10.1007/s11222-022-10097-z](https://doi.org/10.1007/s11222-022-10097-z).

[2] Rothenhäusler, D., Meinshausen, N., Bühlmann, P., & Peters, J. (2021). Anchor regression: Heterogeneous data meet causality. Journal of the Royal Statistical Society Series B: Statistical Methodology, 83(2), 215-246. [doi:10.1111/rssb.12398](https://doi.org/10.1111/rssb.12398).

[3] Kook, L., Baumann, P. F., Dürr, O., Sick, B., & Rügamer, D. (2022). Estimating Conditional Distributions with Neural Networks using R package deeptrafo. arXiv preprint arXiv:2211.13665.
[doi:10.48550/arXiv.2211.13665](https://doi.org/10.48550/arXiv.2211.13665)
