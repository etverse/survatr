## Oracle cross-check of the gcomp survival curve against
## `lmtp::lmtp_tmle(outcome_type = "survival")`. lmtp's TMLE targets the same
## counterfactual survival functional S^a(t); on a confounded DGP, gcomp with
## the correct adjustment set and lmtp both recover the marginal S^a(t) and
## agree up to finite-sample + estimator-form differences. The EIF-based SE is
## NOT comparable to the M-estimation sandwich, so this is a point-estimate
## oracle only. lmtp 1.5.3 returns the survival at the final supplied period, so
## the shared `lmtp_ipw_survival_oracle()` helper loops one fit per horizon and
## reads the S7 `ife@x` estimate.

skip_if_not_installed("lmtp")

test_that("gcomp survival curve agrees with lmtp::lmtp_tmle", {
  skip_on_cran()
  ## Confounded DGP (A depends on L which drives the hazard); gcomp adjusts for
  ## L, so it is correctly specified and recovers the marginal S^a(t).
  dt <- sim_confounded_survival(n = 2000L, K = 5L, seed = 53L, gamma = 0.8)
  fit <- surv_fit(dt, "Y", "A", ~L, "id", "t", time_formula = ~ factor(t))

  tt <- c(3L, 5L)
  s_survatr <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = tt,
    type = "survival"
  )$estimates
  s1 <- s_survatr[get("intervention") == "a1"]$s_hat
  s0 <- s_survatr[get("intervention") == "a0"]$s_hat

  o1 <- lmtp_ipw_survival_oracle(dt, 1, confounder = "L", times = tt)
  o0 <- lmtp_ipw_survival_oracle(dt, 0, confounder = "L", times = tt)
  skip_if(is.null(o1) || is.null(o0), "lmtp_tmle did not fit the toy DGP")

  expect_equal(s1, o1, tolerance = 0.05)
  expect_equal(s0, o0, tolerance = 0.05)
})
