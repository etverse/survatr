test_that("ci_method = 'sandwich' fills se and CI columns on survival curve", {
  dt <- sim_constant_hazard(n = 2000L, K = 8L, h = 0.06, seed = 101L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:8,
    type = "survival",
    ci_method = "sandwich"
  )
  expect_false(any(is.na(res$estimates$se)))
  expect_false(any(is.na(res$estimates$ci_lower)))
  expect_false(any(is.na(res$estimates$ci_upper)))
  expect_true(all(res$estimates$ci_lower <= res$estimates$s_hat))
  expect_true(all(res$estimates$ci_upper >= res$estimates$s_hat))
  expect_true(all(res$estimates$se >= 0))
})

test_that("sandwich CI covers (1 - h)^t on a constant-hazard DGP (single seed)", {
  h <- 0.05
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = h, seed = 113L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = c(1, 5, 10),
    type = "survival",
    ci_method = "sandwich"
  )
  truth <- (1 - h)^c(1, 5, 10)
  covered <- res$estimates$ci_lower <= truth & truth <= res$estimates$ci_upper
  expect_true(all(covered))
})

test_that("sandwich CI coverage ~= nominal across 200 reps", {
  skip_on_cran()
  h <- 0.08
  n <- 1000L
  B <- 200L
  times <- c(1, 5)
  truth <- (1 - h)^times

  covered <- matrix(FALSE, B, length(times))
  for (b in seq_len(B)) {
    dt <- sim_constant_hazard(n = n, K = 6L, h = h, seed = 1000L + b)
    fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
    res <- contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = times,
      type = "survival",
      ci_method = "sandwich"
    )
    covered[b, ] <- res$estimates$ci_lower <= truth &
      truth <= res$estimates$ci_upper
  }
  cov_rate <- colMeans(covered)
  ## Nominal 95%. With B = 200 the binomial SE on a true 0.95 rate is
  ## ~ sqrt(0.95 * 0.05 / 200) ~ 0.015, so the 3-SE band is [0.905, 0.995].
  ## Accept 0.88 as the lower guard to allow for small finite-sample bias.
  expect_true(all(cov_rate >= 0.88))
})
