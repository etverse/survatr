test_that("sandwich log-RR CI is strictly positive and covers 1 under null", {
  dt <- sim_constant_hazard(n = 3000L, K = 8L, h = 0.08, seed = 181L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(1, 4, 8),
    type = "risk_ratio",
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_true(all(res$contrasts$ci_lower > 0))
  expect_true(all(res$contrasts$ci_upper > res$contrasts$ci_lower))
  ## On a no-effect DGP the CI should include 1 at every time.
  expect_true(all(res$contrasts$ci_lower <= 1 & 1 <= res$contrasts$ci_upper))
})
