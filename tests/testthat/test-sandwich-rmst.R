test_that("sandwich SE for rmst is non-negative and CI is finite", {
  dt <- sim_constant_hazard(n = 1500L, K = 8L, h = 0.06, seed = 151L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:8,
    type = "rmst",
    ci_method = "sandwich"
  )
  expect_true(all(res$estimates$se >= 0))
  expect_true(all(is.finite(res$estimates$ci_lower)))
  expect_true(all(is.finite(res$estimates$ci_upper)))
  ## RMST SE is monotone non-decreasing in t because the trapezoidal
  ## quadrature is a cumulative integral of a positive function -- later
  ## times add variance terms, they do not subtract.
  expect_true(all(diff(res$estimates$se) >= -1e-8))
})

test_that("sandwich rmst_difference CI covers 0 on a no-effect DGP", {
  dt <- sim_constant_hazard(n = 3000L, K = 10L, h = 0.06, seed = 163L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:10,
    type = "rmst_difference",
    reference = "a0",
    ci_method = "sandwich"
  )
  ## With no treatment effect in the DGP the RMST difference should be
  ## close to zero; the CI at t = 10 should comfortably cover 0.
  t10 <- res$contrasts[time == 10L]
  expect_true(t10$ci_lower <= 0 && 0 <= t10$ci_upper)
})
