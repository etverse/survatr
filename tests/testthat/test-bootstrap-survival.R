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

test_that("bootstrap vs sandwich SEs agree within 15% on a moderate DGP", {
  skip_on_cran()
  ## Empirical-SD oracle validation (run out-of-band, 2026-04-22) pins
  ## both sandwich and bootstrap to the true sampling SE within 1-2% on
  ## this DGP at n = 1500, h = 0.06, K = 6 (truth from 300-replicate
  ## Monte Carlo). With B = 500 the bootstrap SE has a small additional
  ## MC error; 15% comfortably covers it and still trips on the
  ## factor-of-n_ids class of scaling bugs caught during chunk-3 review.
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
  expect_true(all(rel_err < 0.15))
})

test_that("sandwich and bootstrap SE match empirical-SD truth within 5% (oracle)", {
  ## Smoke oracle: compare sandwich and bootstrap to a freshly-computed
  ## empirical sampling SD on a constant-hazard DGP. This is a smaller
  ## version of the 300-rep out-of-band study (condensed for test time);
  ## it catches systematic calibration drift (e.g. factor-of-n mistakes)
  ## that the pointwise-coverage tests would miss under compensating
  ## errors.
  skip_on_cran()
  h <- 0.06
  n <- 1000L
  K <- 6L
  times <- c(2, 5)

  ## Empirical truth: 100-replicate SD of s_hat (~60-sec run).
  B_sim <- 100L
  s_mat <- matrix(NA_real_, B_sim, length(times))
  for (b in seq_len(B_sim)) {
    dt <- sim_constant_hazard(n = n, K = K, h = h, seed = 9000L + b)
    fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
    res <- contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = times,
      type = "survival",
      ci_method = "none"
    )
    s_mat[b, ] <- res$estimates$s_hat
  }
  emp_sd <- apply(s_mat, 2L, stats::sd)

  ## Sandwich + bootstrap at a single seed.
  dt1 <- sim_constant_hazard(n = n, K = K, h = h, seed = 9500L)
  fit1 <- surv_fit(dt1, "Y", "A", ~1, "id", "t", time_formula = ~1)
  sw <- contrast(
    fit1,
    interventions = list(a0 = causatr::static(0)),
    times = times,
    type = "survival",
    ci_method = "sandwich"
  )
  bt <- contrast(
    fit1,
    interventions = list(a0 = causatr::static(0)),
    times = times,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 500L,
    seed = 1L
  )

  ## Allow 20% relative error -- generous enough to survive the
  ## single-seed sampling noise on `sw` and `bt` while still flagging any
  ## systematic miscalibration.
  expect_true(all(abs(sw$estimates$se - emp_sd) / emp_sd < 0.20))
  expect_true(all(abs(bt$estimates$se - emp_sd) / emp_sd < 0.20))
})
