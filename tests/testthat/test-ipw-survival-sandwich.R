## point-treatment g-computation IPW stacked-EE sandwich. The gold-standard truth for the combined
## (treatment-model + hazard-MSM) variance is the full two-stage nonparametric
## bootstrap, which refits BOTH models per replicate. The stacked sandwich must
## match it; the naive "weights-as-known" SE (the hazard-only block) must not.

test_that("stacked sandwich SE matches the full two-stage bootstrap SE", {
  skip_on_cran()
  dt <- sim_confounded_survival(n = 700L, K = 5L, seed = 7L, gamma = 0.8)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )

  sw <- contrast(
    fit,
    interventions = ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  bs <- contrast(
    fit,
    interventions = ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "bootstrap",
    n_boot = 400L,
    seed = 101L
  )

  ## Per-time risk-difference SE agreement within 15% (bootstrap MC noise at
  ## B = 400). This is the empirical pin on the single-n_ids cross-term scaling.
  expect_equal(sw$contrasts$se, bs$contrasts$se, tolerance = 0.15)

  ## Per-arm survival SE also agrees.
  s_sw <- sw$estimates[get("intervention") == "a1"]$se
  s_bs <- bs$estimates[get("intervention") == "a1"]$se
  expect_equal(s_sw, s_bs, tolerance = 0.15)
})

test_that("the treatment-model correction narrows the SE (block wired)", {
  ## For stabilized weights, accounting for the estimated propensity REDUCES the
  ## variance relative to treating the weights as known (Robins 1999; Hernan et
  ## al. 2000). So the stacked SE must be strictly below the naive hazard-only
  ## SE, and the correction must be non-trivial.
  dt <- sim_confounded_survival(n = 800L, K = 5L, seed = 12L, gamma = 1.0)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  shared <- prepare_sandwich_shared(fit)
  iv <- causatr::static(1)
  times <- 1:5

  stacked <- compute_survival_if_matrix(
    fit,
    iv,
    times,
    shared$prep,
    shared$fit_idx,
    shared$id_vec,
    shared$unique_ids,
    ipw_corr = shared$ipw_corr
  )
  naive <- compute_survival_if_matrix(
    fit,
    iv,
    times,
    shared$prep,
    shared$fit_idx,
    shared$id_vec,
    shared$unique_ids,
    ipw_corr = NULL
  )
  n_ids <- length(shared$unique_ids)
  se_stacked <- sqrt(diag(crossprod(stacked$IF_mat)) / n_ids^2)
  se_naive <- sqrt(diag(crossprod(naive$IF_mat)) / n_ids^2)

  ## Point estimates are identical (the correction touches only the IF).
  expect_equal(stacked$s_hat, naive$s_hat)
  ## Correction is non-trivial and strictly variance-reducing.
  expect_true(all(se_stacked < se_naive))
  expect_true(max(abs(se_stacked - se_naive) / se_naive) > 0.01)
})

test_that("IPW x mgcv::gam baseline hazard: sandwich tracks the bootstrap", {
  skip_on_cran()
  skip_if_not_installed("mgcv")
  ## The gam hazard MSM uses the lpmatrix basis + Vp bread (as in the gcomp gam
  ## path); the IPW treatment-model correction composes on top via prep$X_fit.
  dt <- sim_confounded_survival(n = 900L, K = 8L, seed = 9L, gamma = 0.8)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ s(t, k = 4),
    estimator = "ipw",
    model_fn = mgcv::gam
  )
  sw <- contrast(
    fit,
    interventions = ivs,
    times = 2:8,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  bs <- contrast(
    fit,
    interventions = ivs,
    times = 2:8,
    type = "risk_difference",
    ci_method = "bootstrap",
    n_boot = 250L,
    seed = 9L
  )
  expect_true(all(is.finite(sw$contrasts$se)))
  expect_equal(sw$contrasts$se, bs$contrasts$se, tolerance = 0.15)
})

test_that("IPW sandwich CIs for the risk difference cover the truth", {
  skip_on_cran()
  ## Truth = the marginal RD(t) under the DGP, integrated over L. The default
  ## reference is the first intervention (a1), so the "a0 vs a1" contrast is
  ## risk^0 - risk^1; match that sign. This g-formula truth coincides with the
  ## IPW-MSM probability limit here (the marginal-hazard misspecification bias
  ## of the alpha(t) + A model is < 0.004 at these times, verified out-of-band).
  tt <- c(2L, 4L)
  risk1 <- 1 - true_marginal_survival(1, tt, K = 5L)
  risk0 <- 1 - true_marginal_survival(0, tt, K = 5L)
  true_rd <- risk0 - risk1

  R <- 200L
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  covered <- matrix(FALSE, R, length(tt))
  for (r in seq_len(R)) {
    dt <- sim_confounded_survival(
      n = 800L,
      K = 5L,
      seed = 1000L + r,
      gamma = 0.8
    )
    fit <- surv_fit(
      dt,
      "Y",
      "A",
      ~L,
      "id",
      "t",
      time_formula = ~ factor(t),
      estimator = "ipw"
    )
    res <- contrast(
      fit,
      interventions = ivs,
      times = tt,
      type = "risk_difference",
      ci_method = "sandwich"
    )
    cc <- res$contrasts
    covered[r, ] <- true_rd >= cc$ci_lower & true_rd <= cc$ci_upper
  }
  coverage <- colMeans(covered)
  ## Nominal 95%; allow Monte Carlo slack at R = 200.
  expect_true(all(coverage >= 0.88))
})
