test_that("contrast returns the canonical survatr_result shape", {
  dt <- sim_constant_hazard(n = 2000L, K = 8L, h = 0.06, seed = 5L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:8,
    type = "risk_difference",
    reference = "a0"
  )
  expect_s3_class(res, "survatr_result")
  expect_named(
    res,
    c(
      "estimates",
      "contrasts",
      "time_grid",
      "type",
      "reference",
      "ci_method",
      "call"
    )
  )
  expect_equal(nrow(res$estimates), 16L) ## 2 interventions x 8 times
  expect_equal(nrow(res$contrasts), 8L) ## 1 non-reference x 8 times
  expect_equal(res$time_grid, 1:8)
  expect_equal(res$type, "risk_difference")
  expect_equal(res$reference, "a0")
  expect_equal(res$ci_method, "none")
})

test_that("risk_difference on a DGP with no treatment effect is close to zero", {
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = 0.05, seed = 17L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(1, 5, 10),
    type = "risk_difference",
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, c(0, 0, 0), tolerance = 0.02)
  expect_equal(res$contrasts$contrast, rep("a1 vs a0", 3L))
})

test_that("risk_ratio on a DGP with no treatment effect is close to 1", {
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = 0.05, seed = 19L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(1, 5, 10),
    type = "risk_ratio",
    reference = "a0"
  )
  expect_equal(res$contrasts$estimate, c(1, 1, 1), tolerance = 0.15)
})

test_that("rmst and rmst_difference are curve-shaped", {
  dt <- sim_constant_hazard(n = 2000L, K = 10L, h = 0.05, seed = 23L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)

  rmst_res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:10,
    type = "rmst"
  )
  expect_true("rmst_hat" %in% names(rmst_res$estimates))
  expect_false("s_hat" %in% names(rmst_res$estimates))
  expect_equal(nrow(rmst_res$contrasts), 0L)

  rmst_diff <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:10,
    type = "rmst_difference",
    reference = "a0"
  )
  expect_equal(nrow(rmst_diff$contrasts), 10L)
  expect_equal(rmst_diff$contrasts$estimate, rep(0, 10L), tolerance = 0.1)
})

test_that("all se / ci columns at chunk 2 are NA_real_", {
  dt <- sim_constant_hazard(n = 1000L, K = 5L, h = 0.1, seed = 29L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference"
  )
  expect_true(all(is.na(res$estimates$se)))
  expect_true(all(is.na(res$estimates$ci_lower)))
  expect_true(all(is.na(res$estimates$ci_upper)))
  expect_true(all(is.na(res$contrasts$se)))
})

test_that("contrast is idempotent with respect to fit$pp_data", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 31L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  snap_names <- names(fit$pp_data)
  snap_n <- nrow(fit$pp_data)
  snap_A <- fit$pp_data$A
  contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5
  )
  ## fit$pp_data must not be mutated by the counterfactual intervention.
  expect_identical(names(fit$pp_data), snap_names)
  expect_equal(nrow(fit$pp_data), snap_n)
  expect_equal(fit$pp_data$A, snap_A)
})

test_that("print.survatr_result emits a stable banner", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 37L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference",
    reference = "a0"
  )
  expect_snapshot(print(res))
})
