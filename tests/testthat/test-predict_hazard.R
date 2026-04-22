test_that("predict_hazard_pp returns a vector of length nrow(pp_data)", {
  dt <- sim_constant_hazard(n = 500L, K = 6L, h = 0.08, seed = 1L)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1
  )
  haz <- predict_hazard_pp(fit$model, fit$pp_data)
  expect_type(haz, "double")
  expect_length(haz, nrow(fit$pp_data))
  expect_true(all(haz >= 0 & haz <= 1))
})

test_that("predict_hazard_pp on a constant-hazard DGP recovers the DGP hazard", {
  h_true <- 0.06
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = h_true, seed = 42L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)

  ## Under causatr::static(0) the counterfactual predicted hazard should
  ## be ~ h_true at every row (since the DGP has no treatment / covariate
  ## effect and time_formula = ~ 1).
  pp_cf <- apply_intervention_pp(fit$pp_data, "A", causatr::static(0))
  haz <- predict_hazard_pp(fit$model, pp_cf)
  ## MC-SE on h_hat at n = 5000, K = 10 is ~ sqrt(h(1-h)/(n*K)) ~ 3.4e-3.
  ## Use an absolute tolerance of 0.01 (~3 MC-SE) so the test is robust to
  ## seed-to-seed noise on the DGP.
  expect_lt(abs(mean(haz) - h_true), 0.01)
  expect_equal(stats::sd(haz), 0, tolerance = 1e-8)
})
