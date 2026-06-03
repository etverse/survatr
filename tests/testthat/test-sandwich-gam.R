# GAM (mgcv) hazard models with sandwich variance. The counterfactual design
# matrix must be built on the penalized linear-predictor ("lpmatrix") basis so
# it matches the Bayesian posterior covariance (model$Vp) used as the bread;
# the terms-based model.matrix degrades a smooth term on time to a linear time
# effect and is non-conformable with the bread. These tests guard that path.

test_that("GAM hazard + sandwich runs and returns finite positive SEs", {
  skip_if_not_installed("mgcv")
  dt <- sim_constant_hazard(n = 1000L, K = 8L, h = 0.06, seed = 1L)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ s(t, k = 4),
    model_fn = mgcv::gam
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(2, 4, 6, 8),
    type = "risk_difference",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$contrasts$se)))
  expect_true(all(res$contrasts$se > 0))
})

test_that("GAM sandwich SE matches the analytically-anchored GLM sandwich SE", {
  skip_if_not_installed("mgcv")
  # On a constant-hazard DGP a natural-spline (GLM) and a thin-plate (GAM)
  # baseline fit nearly the same smooth curve, so the standardized
  # risk-difference SEs must coincide. The GLM sandwich is pinned to the
  # closed-form (1 - h)^t variance elsewhere in the sandwich suite, so
  # agreement here transfers that validation to the GAM lpmatrix + Vp path.
  # The comparison is deterministic (no bootstrap Monte-Carlo noise); observed
  # pointwise disagreement is under 2%.
  dt <- sim_constant_hazard(n = 1000L, K = 8L, h = 0.06, seed = 1L)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  tt <- c(2, 4, 6, 8)
  fit_glm <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ splines::ns(t, df = 4),
    model_fn = stats::glm
  )
  fit_gam <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ s(t, k = 4),
    model_fn = mgcv::gam
  )
  se_glm <- contrast(
    fit_glm,
    interventions = ivs,
    times = tt,
    type = "risk_difference",
    ci_method = "sandwich"
  )$contrasts$se
  se_gam <- contrast(
    fit_gam,
    interventions = ivs,
    times = tt,
    type = "risk_difference",
    ci_method = "sandwich"
  )$contrasts$se
  expect_lt(max(abs(se_gam / se_glm - 1)), 0.05)
})
