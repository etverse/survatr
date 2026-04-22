test_that("weights = NULL yields binomial family", {
  dt <- sim_constant_hazard(n = 500L, K = 6L, h = 0.08, seed = 42L)
  fit <- surv_fit(
    data = dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    time_formula = ~1,
    weights = NULL
  )
  expect_equal(fit$family, "binomial")
  expect_equal(fit$model$family$family, "binomial")
})

test_that("non-NULL weights switch the family to quasibinomial", {
  dt <- sim_constant_hazard(n = 500L, K = 6L, h = 0.08, seed = 42L)
  w <- rep(1, nrow(dt))
  fit <- surv_fit(
    data = dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    time_formula = ~1,
    weights = w
  )
  expect_equal(fit$family, "quasibinomial")
  expect_equal(fit$model$family$family, "quasibinomial")
})

test_that("uniform weights ≡ 1 recover the unweighted coefficients", {
  ## Score equations agree between binomial and quasibinomial under w ≡ 1;
  ## the only difference is the dispersion parameter (pinned at 1 vs
  ## estimated). Coefficients must match to machine precision.
  dt <- sim_constant_hazard(n = 1000L, K = 8L, h = 0.06, seed = 7L)
  unw <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1,
    weights = NULL
  )
  w <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1,
    weights = rep(1, nrow(dt))
  )
  expect_equal(
    unname(stats::coef(unw$model)),
    unname(stats::coef(w$model)),
    tolerance = 1e-10
  )
})

test_that("weights validation rejects bad inputs at the surv_fit boundary", {
  dt <- sim_constant_hazard(n = 100L, K = 4L, seed = 1L)
  expect_error(
    surv_fit(dt, "Y", "A", ~1, "id", "t", weights = c(1, 2, 3)),
    class = "survatr_bad_weights"
  )
  expect_error(
    surv_fit(dt, "Y", "A", ~1, "id", "t", weights = rep(-1, nrow(dt))),
    class = "survatr_bad_weights"
  )
})
