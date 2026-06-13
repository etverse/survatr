## Per-cause years of life lost is the integral of the cause-j CIF:
## YLL^(j)(t*) = int_0^{t*} F^(j)(u) du. It reuses the chunk-7 CIF influence
## function mapped through the trapezoidal weights, and satisfies the identity
## sum_j YLL^(j)(t*) = RMTL(t*) (since sum_j F^(j) = 1 - S, partition of unity).
## No `adjustedCurves` oracle is available here; the analytic two-cause CIF
## (`analytic_cr()`) integrated on the grid is the point oracle, and the
## machine-precision identity vs RMTL is the cross-estimand check.

test_that("sum_j YLL^(j)(t*) equals RMTL(t*) to machine precision", {
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  ## Confounded two-cause DGP so the curves are non-trivial.
  dc <- sim_two_cause_constant_hazard(
    n = 2500L,
    K = 20L,
    seed = 401L,
    beta1_A = -0.4,
    beta2_A = 0.3,
    beta1_L = 0.5,
    gamma = 0.6
  )
  fc <- surv_fit(
    dc,
    "event",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "event"
  )
  times <- 1:20

  yll <- contrast(fc, ivs, times, type = "yll", ci_method = "none")$estimates
  rmtl <- contrast(fc, ivs, times, type = "rmtl", ci_method = "none")$estimates

  ## sum over causes per (intervention, time) == all-cause RMTL.
  yll_sum <- yll[,
    list(yll = sum(yll_hat)),
    by = c("intervention", "time")
  ]
  m <- merge(
    yll_sum,
    rmtl[, c("intervention", "time", "rmtl_hat"), with = FALSE],
    by = c("intervention", "time")
  )
  expect_equal(m$yll, m$rmtl_hat, tolerance = 1e-10)
})

test_that("per-cause YLL matches the analytic CIF integral", {
  ## No-effect, no-confounding DGP: the marginal per-cause CIF equals the
  ## closed-form `analytic_cr()`, so the YLL truth is the trapezoidal integral
  ## of that discrete CIF on the grid (same quadrature the estimator uses).
  h1 <- 0.06
  h2 <- 0.04
  dc <- sim_two_cause_constant_hazard(
    n = 6000L,
    K = 25L,
    h1 = h1,
    h2 = h2,
    seed = 409L
  )
  fc <- surv_fit(
    dc,
    "event",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "event"
  )
  times <- 1:25

  yll <- contrast(
    fc,
    list(a0 = causatr::static(0)),
    times,
    type = "yll",
    ci_method = "sandwich"
  )$estimates

  truth_cif <- analytic_cr(times, h1, h2)
  w <- rmst_weights(times) ## trapezoidal weights, F(0) = 0
  truth_yll1 <- as.numeric(w %*% truth_cif$F1)
  truth_yll2 <- as.numeric(w %*% truth_cif$F2)

  est1 <- yll[get("cause") == 1][order(time)]$yll_hat
  est2 <- yll[get("cause") == 2][order(time)]$yll_hat
  ## Cumulative integral over the whole grid; compare at every time.
  expect_equal(est1, truth_yll1, tolerance = 0.05)
  expect_equal(est2, truth_yll2, tolerance = 0.05)
  expect_true(all(yll$se >= 0))
})

test_that("YLL sandwich and bootstrap SE agree", {
  ivs <- list(a0 = causatr::static(0))
  dc <- sim_two_cause_constant_hazard(n = 2000L, K = 15L, seed = 415L)
  fc <- surv_fit(
    dc,
    "event",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "event"
  )
  times <- 1:15

  se_sand <- contrast(
    fc,
    ivs,
    times,
    type = "yll",
    ci_method = "sandwich"
  )$estimates[get("cause") == 1][order(time)]$se
  se_boot <- contrast(
    fc,
    ivs,
    times,
    type = "yll",
    ci_method = "bootstrap",
    n_boot = 120L,
    seed = 3L
  )$estimates[get("cause") == 1][order(time)]$se

  ratio <- se_boot[length(times)] / se_sand[length(times)]
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.5)
})

test_that("YLL on a single-event fit aborts with survatr_yll_needs_cr", {
  dt <- sim_constant_hazard(n = 300L, K = 10L, h = 0.05, seed = 421L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  expect_error(
    contrast(fit, list(a0 = causatr::static(0)), 1:10, type = "yll"),
    class = "survatr_yll_needs_cr"
  )
})
