## Competing-risks CIF: survatr gcomp vs riskRegression::ate g-formula oracle.
## riskRegression::ate fits cause-specific Cox models (CSC) and standardises to
## the observed covariate distribution — the same g-formula estimand as
## survatr's pooled-logistic gcomp. Pooled logistic ≈ Cox for small per-
## interval hazards (h < 0.15), so the two should agree to within 0.03
## absolute CIF units on a correctly-specified confounded DGP at n = 3000.
##
## Oracle: rr_ate_cif_oracle() in helper-cr-oracle.R; requires riskRegression.

test_that("gcomp CIF matches riskRegression::ate on a confounded 2-cause DGP", {
  skip_if_not_installed("riskRegression")

  pp <- sim_two_cause_constant_hazard(
    n = 3000L,
    K = 6L,
    h1 = 0.06,
    h2 = 0.04,
    beta1_A = 0.6,
    beta2_A = -0.4,
    beta1_L = 0.4,
    beta2_L = 0.3,
    gamma = 0.5,
    seed = 7L
  )
  times <- c(2L, 4L, 6L)

  ## survatr gcomp: saturated time formula (correctly specified for the DGP).
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~ factor(t)
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = times,
    type = "cif"
  )

  ## riskRegression oracle: cause-specific Cox + ate() g-formula.
  rr <- rr_ate_cif_oracle(pp, times)

  ## Both estimators target the same marginal standardised CIF; the tolerance
  ## covers the model-class difference (logit vs log-hazard link) + MC error.
  for (j in c(1L, 2L)) {
    for (arm in c("0", "1")) {
      iv_name <- if (arm == "1") "a1" else "a0"
      survatr_cif <- res$estimates[
        get("cause") == j & get("intervention") == iv_name
      ][order(time)]$cif_hat
      rr_cif <- rr[cause == j & as.character(treatment) == arm][order(time)]$cif
      ## Tolerance covers model-class difference (pooled logistic vs Cox) +
      ## MC error; wider than the delicatessen oracle because Cox and logistic
      ## diverge at larger hazards (cause 1, treated arm, ~4% relative).
      expect_equal(
        survatr_cif,
        rr_cif,
        tolerance = 0.06,
        label = sprintf("CIF cause %d arm %s", j, arm)
      )
    }
  }
})
