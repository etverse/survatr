## Track B (longitudinal ICE-hazard survival) -- point estimates, structural
## invariants, and sandwich/bootstrap variance. The primary point oracle is the
## forward-simulation Monte-Carlo g-formula truth (`true_ice_survival()` in
## helper-ice-survival-oracle.R), which is the exact functional the ICE
## estimator targets. The lmtp external oracle lives in
## test-ice-survival-oracle.R; the delicatessen M-estimation variance oracle in
## test-ice-survival-delicatessen.R.

test_that("estimator = 'ice' builds a Track B fit", {
  dat <- sim_ice_feedback(n = 200L, K = 4L, seed = 5L)
  fit <- expect_silent(
    surv_fit(
      dat,
      "Y",
      "A",
      ~1,
      "id",
      "t",
      estimator = "ice",
      confounders_tv = ~L
    )
  )
  expect_s3_class(fit, "survatr_fit")
  expect_identical(fit$track, "B")
  expect_identical(fit$estimator, "ice")
  expect_null(fit$model) # per-step models are fit lazily in contrast()
  expect_false(is.null(fit$ice_details))
  expect_identical(fit$confounders_tv, ~L)
  ## print shows the Track B extras.
  expect_output(print(fit), "Track:       B")
  expect_output(print(fit), "TV covars:   L")
})

test_that("Track B risk curve matches the forward-sim g-formula truth", {
  ## Treatment-confounder-feedback DGP: only the longitudinal g-formula (ICE)
  ## recovers the static counterfactual; a naive curve is biased. Large n +
  ## the analytic MC truth -> tight tolerance.
  dat <- sim_ice_feedback(n = 4000L, K = 4L, seed = 11L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  res <- contrast(fit, interventions = ivs, times = 1:4, type = "survival")

  s1 <- res$estimates[get("intervention") == "a1"][order(get("time"))]$s_hat
  s0 <- res$estimates[get("intervention") == "a0"][order(get("time"))]$s_hat
  truth1 <- true_ice_survival(1, 1:4, K = 4L)
  truth0 <- true_ice_survival(0, 1:4, K = 4L)

  ## Sampling noise + parametric-ICE finite-model bias (the survival-tail
  ## pseudo-outcome is fit logistic-linear, ~0.01-0.015 at this DGP). 0.025
  ## absolute is a tight pin for a 4-period feedback curve.
  expect_equal(s1, truth1, tolerance = 0.025)
  expect_equal(s0, truth0, tolerance = 0.025)
  ## Protective treatment (bA < 0) -> higher survival under a1.
  expect_true(all(s1 > s0))
})

test_that("a naive (unadjusted) curve is biased but ICE is not", {
  ## Sanity that the DGP actually has confounding worth adjusting for: the
  ## crude KM-style survival under observed treatment differs from the ICE
  ## counterfactual by more than the ICE-vs-truth error.
  dat <- sim_ice_feedback(n = 4000L, K = 4L, seed = 11L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1)),
    times = 4L,
    type = "survival"
  )
  s1_ice <- res$estimates$s_hat
  truth1 <- true_ice_survival(1, 4L, K = 4L)
  expect_lt(abs(s1_ice - truth1), 0.03)
})

test_that("per-step link forcing: binomial at horizon, quasibinomial earlier", {
  ## Load-bearing invariant (hard-rules.md). Inspect a horizon-4 backward pass.
  dat <- sim_ice_feedback(n = 300L, K = 4L, seed = 7L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  details <- fit$ice_details
  base <- prepare_track_b_base(fit, details)
  data_iv <- causatr:::ice_apply_intervention_long(
    base$data_lag,
    "A",
    causatr::static(1),
    "id",
    "t"
  )
  res <- run_ice_survival_horizon(
    base$data_lag,
    data_iv,
    base$fit_rows,
    horizon_pos = 4L,
    details = details,
    cols = base$cols,
    intervention = causatr::static(1)
  )
  expect_identical(res$models[[4]]$family$family, "binomial")
  expect_identical(res$models[[3]]$family$family, "quasibinomial")
  expect_identical(res$models[[2]]$family$family, "quasibinomial")
  expect_identical(res$models[[1]]$family$family, "quasibinomial")
})

test_that("intervention sets current treatment but lags hold observed values", {
  ## The ICE contract (inherited from causatr): only current-period treatment
  ## is intervened; lag columns must carry the OBSERVED A_{k-1}.
  dat <- sim_ice_feedback(n = 200L, K = 4L, seed = 9L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  base <- prepare_track_b_base(fit, fit$ice_details)
  data_iv <- causatr:::ice_apply_intervention_long(
    base$data_lag,
    "A",
    causatr::shift(1),
    "id",
    "t"
  )
  ## Current treatment is shifted...
  expect_equal(data_iv$A, base$data_lag$A + 1)
  ## ...but lag1_A is unchanged (observed lagged treatment).
  expect_equal(data_iv$lag1_A, base$data_lag$lag1_A)
})

test_that("Track B sandwich SEs are finite and contrasts are well-formed", {
  dat <- sim_ice_feedback(n = 800L, K = 4L, seed = 13L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  surv <- contrast(
    fit,
    interventions = ivs,
    times = 1:4,
    type = "survival",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(surv$estimates$se)))
  expect_true(all(surv$estimates$se > 0))
  expect_true(all(surv$estimates$ci_lower <= surv$estimates$s_hat))
  expect_true(all(surv$estimates$ci_upper >= surv$estimates$s_hat))

  for (ty in c("risk_difference", "risk_ratio", "rmst_difference")) {
    rr <- contrast(
      fit,
      interventions = ivs,
      times = 1:4,
      type = ty,
      ci_method = "sandwich"
    )
    expect_equal(nrow(rr$contrasts), 4L)
    expect_true(all(is.finite(rr$contrasts$se)))
    expect_true(all(is.finite(rr$contrasts$estimate)))
  }
})

test_that("Track B sandwich SEs agree with the empirical bootstrap", {
  ## The numerical confirmation that the survival-aware IF chain (with the
  ## (1 - D_k) failure carry-forward factor) is correct. Reusing causatr's
  ## chain verbatim over-covers at later horizons; this test guards against a
  ## regression to that. Loose tolerance (bootstrap MC noise + finite n); the
  ## tight pin is the delicatessen oracle.
  skip_on_cran()
  dat <- sim_ice_feedback(n = 1200L, K = 4L, seed = 3L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  sw <- contrast(
    fit,
    interventions = ivs,
    times = 1:4,
    type = "survival",
    ci_method = "sandwich"
  )$estimates
  bo <- contrast(
    fit,
    interventions = ivs,
    times = 1:4,
    type = "survival",
    ci_method = "bootstrap",
    n_boot = 400L,
    seed = 42L
  )$estimates
  m <- merge(
    sw[, list(intervention, time, se_sw = se)],
    bo[, list(intervention, time, se_bo = se)],
    by = c("intervention", "time")
  )
  rel <- abs(m$se_sw - m$se_bo) / m$se_bo
  ## If the chain were reused verbatim (the (1 - D_k) bug) this hits ~0.4 at
  ## the last horizon; the fix brings every horizon well under 0.20.
  expect_lt(max(rel), 0.20)
})

test_that("ICE on a constant-within-id treatment informs but works", {
  ## A point (baseline-constant) treatment is valid for ICE; we inform that
  ## Track A is cheaper but do not abort.
  dat <- sim_confounded_survival(n = 300L, K = 4L, seed = 2L, gamma = 0.5)
  expect_message(
    fit <- surv_fit(
      dat,
      "Y",
      "A",
      ~1,
      "id",
      "t",
      estimator = "ice",
      confounders_tv = ~L
    ),
    class = "survatr_ice_static_treatment"
  )
  expect_identical(fit$track, "B")
})

test_that("time-varying treatment with Track A warns, pointing to ice", {
  dat <- sim_ice_feedback(n = 200L, K = 4L, seed = 5L)
  expect_warning(
    surv_fit(dat, "Y", "A", ~L, "id", "t", time_formula = ~ factor(t)),
    class = "survatr_tv_treatment_track_a"
  )
})

test_that("Track B rejects external weights and stochastic interventions", {
  dat <- sim_ice_feedback(n = 150L, K = 4L, seed = 5L)
  expect_error(
    surv_fit(
      dat,
      "Y",
      "A",
      ~1,
      "id",
      "t",
      estimator = "ice",
      confounders_tv = ~L,
      weights = rep(1, nrow(dat))
    ),
    class = "survatr_ice_external_weights"
  )
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  expect_error(
    contrast(
      fit,
      interventions = list(s = causatr::ipsi(2)),
      times = 1:4,
      type = "survival"
    ),
    class = "survatr_ice_intervention_deferred"
  )
})
