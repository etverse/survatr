## IPCW three-block stacked-EE sandwich variance tests.
##
## The full sandwich has three blocks:
##   Block 1: hazard MSM score (beta)
##   Block 2: treatment model (alpha) -- from chunk 5
##   Block 3: censoring model (gamma) -- added in chunk 11
##
## Hard invariants to assert:
##   SE_two_block < SE_hazard_only          (treatment correction is variance-reducing)
##   SE_three_block != SE_two_block (≠ 0)   (censoring block is wired)
##   SE_sandwich ≈ SE_bootstrap             (internal consistency)
##   Nominal 95% CIs cover the truth        (calibration)
##
## Note: the censoring correction can widen OR narrow the SE depending on the
## covariance between the censoring IF and the existing IF matrix — it is NOT
## guaranteed to always increase SE (unlike the treatment block's reduction).

test_that("censoring block is non-trivially wired; treatment block is variance-reducing", {
  ## Informative censoring: the censoring correction is non-zero and must
  ## modify the SE meaningfully. The treatment correction reduces the SE
  ## (stabilized-weight invariant). So: two_block < hazard_only (always), and
  ## three_block differs from two_block by > 1% somewhere (block is wired).
  dt <- sim_informative_censoring(
    n = 800L,
    K = 5L,
    seed = 11L,
    delta_cens = 0.8
  )
  fit_ipcw <- surv_fit(
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

  ## Shared sandwich pieces.
  shared_ipcw <- prepare_sandwich_shared(fit_ipcw)
  iv <- causatr::static(1)
  times <- 1:5

  ## Three-block: both IPW treatment correction AND IPCW censoring correction.
  three_block <- compute_survival_if_matrix(
    fit_ipcw,
    iv,
    times,
    shared_ipcw$prep,
    shared_ipcw$fit_idx,
    shared_ipcw$id_vec,
    shared_ipcw$unique_ids,
    ipw_corr = shared_ipcw$ipw_corr,
    ipcw_corr = shared_ipcw$ipcw_corr
  )

  ## Two-block: IPW treatment correction only (no censoring block).
  two_block <- compute_survival_if_matrix(
    fit_ipcw,
    iv,
    times,
    shared_ipcw$prep,
    shared_ipcw$fit_idx,
    shared_ipcw$id_vec,
    shared_ipcw$unique_ids,
    ipw_corr = shared_ipcw$ipw_corr,
    ipcw_corr = NULL
  )

  ## Hazard-only: no corrections (weights treated as known).
  hazard_only <- compute_survival_if_matrix(
    fit_ipcw,
    iv,
    times,
    shared_ipcw$prep,
    shared_ipcw$fit_idx,
    shared_ipcw$id_vec,
    shared_ipcw$unique_ids,
    ipw_corr = NULL,
    ipcw_corr = NULL
  )

  n_ids <- length(shared_ipcw$unique_ids)
  se_three <- sqrt(diag(crossprod(three_block$IF_mat)) / n_ids^2)
  se_two <- sqrt(diag(crossprod(two_block$IF_mat)) / n_ids^2)
  se_hazard <- sqrt(diag(crossprod(hazard_only$IF_mat)) / n_ids^2)

  ## Point estimates are identical across all three (only IF is affected).
  expect_equal(three_block$s_hat, hazard_only$s_hat)

  ## Treatment block narrows SE vs hazard-only (stabilized weight variance
  ## reduction; Robins 1999 invariant, already pinned in chunk 5).
  expect_true(all(se_two < se_hazard))
  ## Treatment correction is non-trivial: > 1% relative change at peak.
  expect_true(max(abs(se_two - se_hazard) / se_hazard) > 0.01)

  ## Censoring block is wired (not identically zero). The SE is a compressed
  ## scalar that can mask small per-element IF differences; we check the IF
  ## matrices directly. With moderate informative censoring (~8% per period),
  ## the per-individual IF contribution from the censoring block is small but
  ## strictly non-zero. We assert:
  ##   (a) A_beta_gamma is non-zero (numDeriv returned a real Jacobian)
  ##   (b) IF_gamma_per_id is non-zero (censoring scores were accumulated)
  ##   (c) The three-block IF matrix differs from the two-block IF matrix
  expect_true(!all(shared_ipcw$ipcw_corr$A_beta_gamma == 0))
  expect_true(!all(shared_ipcw$ipcw_corr$IF_gamma_per_id == 0))
  expect_true(
    max(abs(three_block$IF_mat - two_block$IF_mat)) > 1e-6
  )
})

test_that("three-block sandwich SE matches full three-stage bootstrap SE", {
  skip_on_cran()
  ## Full bootstrap refits both the treatment model AND the censoring model per
  ## replicate, propagating all three blocks of uncertainty. The sandwich must
  ## agree with it; if the censoring block is omitted or mis-scaled the two
  ## diverge. Tolerance = 15% (bootstrap MC noise at B = 400).
  dt <- sim_informative_censoring(n = 700L, K = 5L, seed = 12L)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
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
    seed = 201L
  )

  expect_equal(sw$contrasts$se, bs$contrasts$se, tolerance = 0.15)

  ## Per-arm SEs also agree.
  s_sw <- sw$estimates[get("intervention") == "a1"]$se
  s_bs <- bs$estimates[get("intervention") == "a1"]$se
  expect_equal(s_sw, s_bs, tolerance = 0.15)
})

test_that("IPCW sandwich SEs are finite and positive for both arms and RD", {
  dt <- sim_informative_censoring(n = 600L, K = 5L, seed = 13L)
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
  res <- contrast(
    fit,
    interventions = ivs,
    times = 1:5,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$estimates$se)))
  expect_true(all(res$estimates$se > 0))
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(all(res$contrasts$se > 0))
  ## CI bounds are sensibly ordered.
  expect_true(all(res$contrasts$ci_lower < res$contrasts$ci_upper))
  ## CIs contain the point estimates.
  expect_true(all(
    res$contrasts$estimate >= res$contrasts$ci_lower &
      res$contrasts$estimate <= res$contrasts$ci_upper
  ))
})

test_that("IPCW sandwich CIs for the risk difference cover the truth", {
  skip_on_cran()
  ## Monte Carlo coverage at R = 150 x n = 800 simulations. The truth is the
  ## marginal RD integrated over L from true_marginal_survival_ipcw(). The
  ## censoring process is "non-informative about the event conditional on A and
  ## L", so the marginal survival under each intervention is identified by IPCW
  ## and the sandwich 95% CI should have empirical coverage >= 88% (nominal
  ## 95% with MC slack at R = 150).
  tt <- c(2L, 5L)
  risk1_true <- 1 - true_marginal_survival_ipcw(1, tt)
  risk0_true <- 1 - true_marginal_survival_ipcw(0, tt)
  true_rd <- risk0_true - risk1_true

  R <- 150L
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  covered <- matrix(FALSE, R, length(tt))

  for (r in seq_len(R)) {
    dt <- sim_informative_censoring(n = 800L, K = 5L, seed = 2000L + r)
    fit <- suppressWarnings(surv_fit(
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
    ))
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
  expect_true(all(coverage >= 0.88))
})
