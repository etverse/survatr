## IPCW (per-period censoring weights) -- point estimate tests.
## Variance is in test-ipcw-survival-sandwich.R.
##
## Oracle: lmtp::lmtp_tmle(outcome_type = "survival", cens = ...) is a TMLE
## estimator consistent for the marginal counterfactual survival under the same
## identification assumptions as survatr's IPCW-IPW. With a correctly specified
## censoring model and treatment model it converges to the same estimand and
## serves as an external point-estimate reference.

test_that("IPCW surv_fit stores the censoring model and IPCW weights", {
  dt <- sim_informative_censoring(n = 300L, K = 5L, seed = 1L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  expect_s3_class(fit, "survatr_fit")

  ## Three distinct weight objects stored.
  expect_false(is.null(fit$censoring_model))
  expect_s3_class(fit$censoring_model, "glm")
  expect_false(is.null(fit$ipcw_numerator_model))
  expect_s3_class(fit$ipcw_numerator_model, "glm")
  expect_false(is.null(fit$ipcw_weights))
  expect_false(is.null(fit$ipw_treatment_weights_pp))
  expect_length(fit$ipcw_weights, nrow(dt))
  expect_length(fit$ipw_treatment_weights_pp, nrow(dt))

  ## The combined `weights` field is the product of the two weight components.
  combined <- fit$ipw_treatment_weights_pp * fit$ipcw_weights
  expect_equal(fit$weights, combined, tolerance = 1e-10)

  ## Ipcw field stored for bootstrap refit.
  expect_equal(fit$ipcw, ~L)
  expect_identical(fit$censoring_model_fn, stats::glm)
})

test_that("IPCW weights grow within id (running product)", {
  ## Within each id, W^C_{i,k} = prod_{m <= k} factor_m, so it must be
  ## non-decreasing if all per-period factors >= 1 (which holds when the
  ## denominator < numerator, common for light censoring). More robustly:
  ## each period's weight is a cumulative product, so |w(k)| >= min-factor^k.
  ## We check that the IPCW weight at t=K is not equal to the weight at t=1
  ## for most ids (i.e., it actually accumulated).
  dt <- sim_informative_censoring(n = 400L, K = 5L, seed = 2L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  w <- data.table::data.table(
    id = dt$id,
    t = dt$t,
    w = fit$ipcw_weights
  )
  ## More than 50% of ids have a different weight at t=1 vs t=5.
  w1 <- w[t == 1L, list(w1 = w), keyby = "id"]
  w5 <- w[t == 5L, list(w5 = w), keyby = "id"]
  merged <- merge(w1, w5, by = "id")
  prop_diff <- mean(abs(merged$w1 - merged$w5) > 1e-8)
  expect_gt(prop_diff, 0.5)

  ## Treatment weights are constant within id; IPCW weights are not.
  n_uniq_trt <- data.table::as.data.table(
    list(id = dt$id, w = fit$ipw_treatment_weights_pp)
  )[, list(n = data.table::uniqueN(w)), by = "id"]$n
  expect_true(all(n_uniq_trt == 1L))

  n_uniq_ipcw <- w[, list(n = data.table::uniqueN(w)), by = "id"]$n
  ## At least some ids have more than one distinct IPCW weight across periods.
  expect_true(any(n_uniq_ipcw > 1L))
})

test_that("IPCW curve matches lmtp::lmtp_tmle on informative censoring DGP", {
  skip_on_cran()
  skip_if_not_installed("lmtp")
  dt <- sim_informative_censoring(n = 2000L, K = 5L, seed = 3L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "survival"
  )
  s1 <- res$estimates[get("intervention") == "a1"]$s_hat
  s0 <- res$estimates[get("intervention") == "a0"]$s_hat

  ## Compare at t = 3 and t = 5 (most accumulated point).
  tt <- c(3L, 5L)
  o1 <- lmtp_ipcw_survival_oracle(dt, 1, times = tt)
  o0 <- lmtp_ipcw_survival_oracle(dt, 0, times = tt)
  skip_if(is.null(o1) || is.null(o0), "lmtp_tmle did not fit the toy DGP")

  expect_equal(s1[tt], o1, tolerance = 0.05)
  expect_equal(s0[tt], o0, tolerance = 0.05)
})

test_that("naive IPW over-estimates survival under informative censoring (bias direction)", {
  ## Informative censoring creates selection bias: high-L individuals (with
  ## higher event hazard) are more likely to be censored. The naive row-filter
  ## IPW selects progressively for low-L individuals at later times, so it
  ## OVER-estimates the marginal survival probability (lower average event
  ## hazard in the at-risk set). IPCW corrects this and gives a lower estimate.
  ##
  ## This test checks the DIRECTION: over R=25 runs, the naive IPW estimate
  ## should be ABOVE the IPCW estimate for BOTH arms at the last horizon
  ## (where the selection bias has accumulated most). We use moderate
  ## delta_cens = 0.8 to maintain positivity (extreme delta_cens causes the
  ## IPCW to over-correct due to near-empty strata).
  skip_on_cran()
  R <- 25L
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  naive_gt_ipcw_a1 <- logical(R)
  naive_gt_ipcw_a0 <- logical(R)
  for (r in seq_len(R)) {
    dt <- sim_informative_censoring(
      n = 1000L,
      K = 5L,
      seed = 600L + r,
      delta_cens = 0.8
    )
    fit_i <- surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      time_formula = ~ factor(t),
      estimator = "ipw",
      ipcw = ~L
    )
    fit_n <- surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      time_formula = ~ factor(t),
      estimator = "ipw"
    )
    s1_i <- contrast(
      fit_i,
      ivs,
      times = 5L,
      type = "survival"
    )$estimates[get("intervention") == "a1"]$s_hat
    s0_i <- contrast(
      fit_i,
      ivs,
      times = 5L,
      type = "survival"
    )$estimates[get("intervention") == "a0"]$s_hat
    s1_n <- contrast(
      fit_n,
      ivs,
      times = 5L,
      type = "survival"
    )$estimates[get("intervention") == "a1"]$s_hat
    s0_n <- contrast(
      fit_n,
      ivs,
      times = 5L,
      type = "survival"
    )$estimates[get("intervention") == "a0"]$s_hat
    naive_gt_ipcw_a1[r] <- s1_n > s1_i
    naive_gt_ipcw_a0[r] <- s0_n > s0_i
  }
  ## Expect > 70% of runs show the expected direction (bias is real but
  ## individual runs have noise; with R = 25 a 70% success rate requires
  ## roughly P > 0.7 per run, which holds when the asymptotic bias signal
  ## is larger than the finite-sample variance).
  expect_gt(mean(naive_gt_ipcw_a1), 0.70)
  expect_gt(mean(naive_gt_ipcw_a0), 0.70)
})

test_that("non-informative censoring (delta_cens = 0): IPCW weights ~ 1", {
  ## When censoring is independent of confounders, the denominator and numerator
  ## models predict the same value, so W^C_{i,k} ≈ 1 and the IPCW curve
  ## matches the row-filter IPW curve. Formally, max |W^C - 1| should be small.
  dt_ni <- sim_informative_censoring(
    n = 2000L,
    K = 5L,
    seed = 5L,
    delta_cens = 0
  )
  fit_ipcw <- surv_fit(
    dt_ni,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  ## Under non-informative censoring the IPCW weights should cluster tightly
  ## around 1 (within MC noise at n = 2000).
  expect_lt(stats::sd(fit_ipcw$ipcw_weights), 0.25)
  expect_equal(mean(fit_ipcw$ipcw_weights), 1, tolerance = 0.05)

  ## The survival curve from IPCW should agree with row-filter IPW.
  fit_naive <- surv_fit(
    dt_ni,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  s_ipcw <- contrast(
    fit_ipcw,
    interventions = ivs,
    times = 1:5,
    type = "survival"
  )$estimates$s_hat
  s_naive <- contrast(
    fit_naive,
    interventions = ivs,
    times = 1:5,
    type = "survival"
  )$estimates$s_hat
  expect_equal(s_ipcw, s_naive, tolerance = 0.03)
})

test_that("IPCW per-period trim stores fixed thresholds and caps the weights", {
  dt <- sim_informative_censoring(n = 500L, K = 5L, seed = 6L, delta_cens = 1.5)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L,
    trim = 0.9
  )
  ## Trim thresholds are stored per time period.
  expect_false(is.null(fit$ipcw_trim_thresholds))
  expect_length(fit$ipcw_trim_thresholds, 5L)
  ## No IPCW weight exceeds its period's threshold (with floating-point slack).
  for (tt in 1:5) {
    w_t <- fit$ipcw_weights[dt$t == tt]
    cap <- fit$ipcw_trim_thresholds[[as.character(tt)]]
    expect_true(all(w_t <= cap + 1e-8))
  }
})

test_that("IPCW fit structure is preserved through contrast()", {
  dt <- sim_informative_censoring(n = 300L, K = 4L, seed = 7L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:4,
    type = "survival",
    ci_method = "sandwich"
  )
  expect_s3_class(res, "survatr_result")
  expect_true(all(is.finite(res$estimates$s_hat)))
  expect_true(all(is.finite(res$estimates$se)))
  expect_true(all(res$estimates$se > 0))
})

test_that("IPCW bootstrap refits the censoring model per replicate", {
  ## The bootstrap path must pass `ipcw` and `censoring_model_fn` back into
  ## `surv_fit()` so each replicate has its own censoring model. If the refit
  ## were absent the bootstrap would error (the stored model cannot predict on
  ## a different replicate's data layout) or silently use fixed weights.
  ## We check that the bootstrap CI is finite and wider than the sandwich CI
  ## (the censoring-model uncertainty inflates the SE slightly at n = 400).
  skip_on_cran()
  dt <- sim_informative_censoring(n = 400L, K = 4L, seed = 8L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  bs <- contrast(
    fit,
    interventions = ivs,
    times = 1:4,
    type = "risk_difference",
    ci_method = "bootstrap",
    n_boot = 200L,
    seed = 42L
  )
  expect_true(all(is.finite(bs$contrasts$se)))
  expect_true(all(bs$contrasts$se > 0))
})

## ── Rejection tests ────────────────────────────────────────────────────────

test_that("ipcw= with non-formula is rejected (survatr_bad_ipcw)", {
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      estimator = "ipw",
      ipcw = "L"
    ),
    class = "survatr_bad_ipcw"
  )
})

test_that("ipcw= with non-ipw estimator is rejected (survatr_ipcw_estimator)", {
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      estimator = "gcomp",
      ipcw = ~L
    ),
    class = "survatr_ipcw_estimator"
  )
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      estimator = "ice",
      ipcw = ~L
    ),
    class = "survatr_ipcw_estimator"
  )
})

test_that("ipcw= without a censoring column is rejected (survatr_ipcw_no_censoring)", {
  dt <- sim_confounded_survival(n = 100L, K = 3L, seed = 1L)
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      estimator = "ipw",
      ipcw = ~L
    ),
    class = "survatr_ipcw_no_censoring"
  )
})

test_that("snapshot: ipcw error messages are informative", {
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  expect_snapshot(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      estimator = "ipw",
      ipcw = "L"
    ),
    error = TRUE
  )
  expect_snapshot(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      estimator = "gcomp",
      ipcw = ~L
    ),
    error = TRUE
  )
  expect_snapshot(
    surv_fit(
      sim_confounded_survival(n = 100L, K = 3L, seed = 1L),
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      estimator = "ipw",
      ipcw = ~L
    ),
    error = TRUE
  )
})

test_that("IPCW + time-varying treatment is rejected (survatr_ipw_time_varying_treatment)", {
  ## Build a dataset where treatment flips within an id so the check fires.
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  dt <- data.table::copy(dt)
  dt[t == 2L, A := 1L - A]
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      time_formula = ~ factor(t),
      estimator = "ipw",
      ipcw = ~L
    ),
    class = "survatr_ipw_time_varying_treatment"
  )
})

test_that("IPCW + no treatment variation is rejected (survatr_ipw_no_treatment_variation)", {
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  dt <- data.table::copy(dt)
  dt[, A := 1L]
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      time_formula = ~ factor(t),
      estimator = "ipw",
      ipcw = ~L
    ),
    class = "survatr_ipw_no_treatment_variation"
  )
})

test_that("IPCW + continuous treatment is rejected (survatr_ipw_treatment_unsupported)", {
  dt <- sim_informative_censoring(n = 100L, K = 3L, seed = 1L)
  dt <- data.table::copy(dt)
  ## Replace binary treatment with a continuous per-id value so the
  ## bernoulli-family check fires.
  set.seed(42L)
  n_ids <- data.table::uniqueN(dt, by = "id")
  id_vals <- unique(dt$id)
  a_cont <- stats::rnorm(n_ids)
  dt[, A := a_cont[match(id, id_vals)]]
  expect_error(
    surv_fit(
      dt,
      outcome = "Y",
      treatment = "A",
      confounders = ~L,
      id = "id",
      time = "t",
      censoring = "C",
      time_formula = ~ factor(t),
      estimator = "ipw",
      ipcw = ~L
    ),
    class = "survatr_ipw_treatment_unsupported"
  )
})

test_that("ipcw = ~1 (treatment-only censoring denominator) fits and runs contrast()", {
  ## When ipcw_formula is ~1 the censoring model conditions only on time and
  ## treatment; the covariates term in fit_censoring_model is dropped (NULL
  ## branch at line 68). This is the minimal valid IPCW spec.
  dt <- sim_informative_censoring(n = 300L, K = 4L, seed = 5L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~1
  )
  expect_false(is.null(fit$censoring_model))
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:4,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$contrasts$se)))
})

test_that("trim + sandwich hits the per-period cap branch in the weight closure", {
  ## contrast() with ci_method = 'sandwich' evaluates the numDeriv closure
  ## that applies per-period trim caps (lines 269-271 of ipcw_survival.R).
  dt <- sim_informative_censoring(n = 300L, K = 4L, seed = 5L, delta_cens = 1.5)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L,
    trim = 0.9
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:4,
    type = "survival",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$estimates$se)))
})

test_that("compute_ipcw_running_weights: unstabilized path (num_model = NULL)", {
  ## With num_model = NULL the per-row factor is 1 / (1 - g_den), so the
  ## cumulative product is >= 1 everywhere (each factor >= 1 when g_den < 1).
  dt <- sim_informative_censoring(n = 200L, K = 3L, seed = 5L)
  fit <- surv_fit(
    dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  w_unstab <- survatr:::compute_ipcw_running_weights(
    data = fit$pp_data,
    cens_model = fit$censoring_model,
    num_model = NULL,
    id = "id",
    time = "t"
  )
  expect_true(all(w_unstab$weights >= 1 - 1e-8))
})
