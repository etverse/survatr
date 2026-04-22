## Closed-form oracle for the pooled-logistic family switch.
##
## DGP: constant discrete-time hazard `h` on the at-risk population, no
## covariate or treatment effect, time_formula = ~ 1. The MLE of the logit
## intercept converges to `qlogis(h)` as n -> infty. With n = 5000 the MC
## standard error on beta_0 is roughly sqrt(1 / (n * K * h * (1 - h))) for
## the full at-risk subset; at h = 0.05, K = 10, that's ~0.06 — tolerance 0.1
## gives ~1.6 MC-SE cushion, safe against the binary rbinom noise.

test_that("binomial intercept recovers qlogis(h) on constant-hazard DGP", {
  h <- 0.05
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = h, seed = 123L)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1,
    weights = NULL
  )
  co <- stats::coef(fit$model)
  expect_equal(unname(co[["(Intercept)"]]), stats::qlogis(h), tolerance = 0.1)
  ## No treatment effect in the DGP — beta_A should be close to zero.
  expect_lt(abs(unname(co[["A"]])), 0.2)
})

test_that("quasibinomial intercept recovers qlogis(h) with uniform weights", {
  h <- 0.08
  dt <- sim_constant_hazard(n = 5000L, K = 8L, h = h, seed = 321L)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1,
    weights = rep(1, nrow(dt))
  )
  co <- stats::coef(fit$model)
  expect_equal(unname(co[["(Intercept)"]]), stats::qlogis(h), tolerance = 0.1)
  expect_lt(abs(unname(co[["A"]])), 0.2)
})
