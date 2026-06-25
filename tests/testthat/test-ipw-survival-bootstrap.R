## point-treatment g-computation IPW bootstrap. Each replicate refits BOTH the treatment model and the
## weighted hazard MSM (the weights are re-estimated, not carried), so the
## bootstrap is a genuine two-stage resample. Reproducibility + populated CIs +
## agreement with the sandwich are the wiring checks (the SE-vs-sandwich pin
## lives in test-ipw-survival-sandwich.R).

test_that("IPW bootstrap populates CIs and brackets the point estimate", {
  skip_on_cran()
  dt <- sim_confounded_survival(n = 500L, K = 5L, seed = 21L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  res <- contrast(
    fit,
    interventions = ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "bootstrap",
    n_boot = 200L,
    seed = 5L
  )

  cc <- res$contrasts
  expect_false(anyNA(cc$se))
  expect_false(anyNA(cc$ci_lower))
  expect_true(all(cc$ci_lower <= cc$estimate & cc$estimate <= cc$ci_upper))
  expect_true(all(cc$se > 0))
})

test_that("IPW bootstrap is reproducible under a fixed seed", {
  skip_on_cran()
  dt <- sim_confounded_survival(n = 400L, K = 4L, seed = 22L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  run <- function() {
    contrast(
      fit,
      interventions = ivs,
      times = 1:4,
      type = "risk_difference",
      ci_method = "bootstrap",
      n_boot = 150L,
      seed = 314L
    )$contrasts
  }
  a <- run()
  b <- run()
  expect_equal(a$se, b$se)
  expect_equal(a$ci_lower, b$ci_lower)
  expect_equal(a$ci_upper, b$ci_upper)
})

test_that("IPW bootstrap survival curve SE is finite and ordered", {
  skip_on_cran()
  dt <- sim_confounded_survival(n = 500L, K = 5L, seed = 23L, gamma = 0.6)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1)),
    times = 1:5,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 200L,
    seed = 7L
  )
  est <- res$estimates
  expect_false(anyNA(est$se))
  expect_true(all(est$se > 0))
  expect_true(all(est$ci_lower <= est$s_hat & est$s_hat <= est$ci_upper))
})
