test_that("sandwich fills contrasts$se / ci for risk_difference", {
  dt <- sim_constant_hazard(n = 2000L, K = 6L, h = 0.06, seed = 131L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:6,
    type = "risk_difference",
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_false(any(is.na(res$contrasts$se)))
  expect_false(any(is.na(res$contrasts$ci_lower)))
  expect_false(any(is.na(res$contrasts$ci_upper)))
  expect_true(all(res$contrasts$ci_lower <= res$contrasts$estimate))
  expect_true(all(res$contrasts$ci_upper >= res$contrasts$estimate))
})

test_that("sandwich RD CI covers 0 on a no-effect DGP at nominal 95%", {
  skip_on_cran()
  B <- 200L
  times <- c(1, 5, 10)
  covered <- matrix(FALSE, B, length(times))
  for (b in seq_len(B)) {
    dt <- sim_constant_hazard(n = 1000L, K = 10L, h = 0.06, seed = 2000L + b)
    fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
    res <- contrast(
      fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = times,
      type = "risk_difference",
      reference = "a0",
      ci_method = "sandwich"
    )
    covered[b, ] <- res$contrasts$ci_lower <= 0 & 0 <= res$contrasts$ci_upper
  }
  cov_rate <- colMeans(covered)
  expect_true(all(cov_rate >= 0.88))
})
