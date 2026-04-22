test_that("ci_method = 'bootstrap' fills se and CI columns", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 201L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:5,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 100L,
    seed = 1L
  )
  expect_false(any(is.na(res$estimates$se)))
  expect_false(any(is.na(res$estimates$ci_lower)))
  expect_false(any(is.na(res$estimates$ci_upper)))
  expect_true(all(res$estimates$ci_lower <= res$estimates$s_hat))
  expect_true(all(res$estimates$ci_upper >= res$estimates$s_hat))
  expect_equal(res$ci_method, "bootstrap")
})

test_that("bootstrap percentile CI covers (1 - h)^t on a single-seed DGP", {
  h <- 0.08
  dt <- sim_constant_hazard(n = 2000L, K = 6L, h = h, seed = 211L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = c(1, 3, 6),
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 300L,
    seed = 2L
  )
  truth <- (1 - h)^c(1, 3, 6)
  expect_true(all(
    res$estimates$ci_lower <= truth & truth <= res$estimates$ci_upper
  ))
})

test_that("bootstrap is reproducible with a fixed seed", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 223L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  r1 <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:5,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 50L,
    seed = 42L
  )
  r2 <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:5,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 50L,
    seed = 42L
  )
  expect_equal(r1$estimates$se, r2$estimates$se)
  expect_equal(r1$estimates$ci_lower, r2$estimates$ci_lower)
})

test_that("bootstrap supports risk_difference with populated contrasts CIs", {
  dt <- sim_constant_hazard(n = 800L, K = 5L, h = 0.08, seed = 229L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference",
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 150L,
    seed = 3L
  )
  expect_false(any(is.na(res$contrasts$se)))
  expect_true(all(res$contrasts$ci_lower <= res$contrasts$estimate))
  expect_true(all(res$contrasts$ci_upper >= res$contrasts$estimate))
})

test_that("bootstrap supports risk_ratio with strictly-positive percentile CIs", {
  dt <- sim_constant_hazard(n = 1000L, K = 5L, h = 0.1, seed = 233L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(1, 3, 5),
    type = "risk_ratio",
    reference = "a0",
    ci_method = "bootstrap",
    n_boot = 150L,
    seed = 5L
  )
  expect_true(all(res$contrasts$ci_lower > 0))
  expect_true(all(res$contrasts$ci_upper > res$contrasts$ci_lower))
})

test_that("bootstrap vs sandwich SEs agree within 30% on a moderate DGP", {
  skip_on_cran()
  dt <- sim_constant_hazard(n = 1500L, K = 6L, h = 0.08, seed = 239L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  sw <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = c(1, 3, 6),
    type = "survival",
    ci_method = "sandwich"
  )
  bt <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = c(1, 3, 6),
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 500L,
    seed = 7L
  )
  rel_err <- abs(sw$estimates$se - bt$estimates$se) / sw$estimates$se
  expect_true(all(rel_err < 0.30))
})
