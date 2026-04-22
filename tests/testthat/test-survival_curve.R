## Closed-form oracle: constant discrete-time hazard h implies
##   S(t) = (1 - h)^t
## on the integer time grid. The `compute_survival_curve()` result should
## converge to this curve as n -> infty; at n = 5000, tolerance 0.01 leaves
## roughly 2 MC-SE of headroom.

test_that("survival curve converges to (1 - h)^t on a constant-hazard DGP", {
  h <- 0.05
  dt <- sim_constant_hazard(n = 5000L, K = 12L, h = h, seed = 11L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)

  pp_cf <- apply_intervention_pp(fit$pp_data, "A", causatr::static(0))
  haz <- predict_hazard_pp(fit$model, pp_cf)

  times <- c(1, 5, 10, 12)
  curve <- compute_survival_curve(
    pp_data = pp_cf,
    hazards = haz,
    id = "id",
    time = "t",
    times = times,
    intervention_name = "a0"
  )

  expect_equal(curve$time, times)
  ## Absolute tolerance: at n = 5000, S_bar has MC-SE ~ sqrt(S(1-S)/n) which
  ## peaks near 0.007; allow 0.03 to cover ~4 MC-SE plus small GLM-fit bias
  ## near the upper end of the time grid.
  expect_lt(max(abs(curve$s_hat - (1 - h)^times)), 0.03)
  expect_lt(max(abs(curve$risk_hat - (1 - (1 - h)^times))), 0.03)
  expect_equal(curve$intervention, rep("a0", length(times)))
  ## SE columns are placeholders at chunk 2.
  expect_true(all(is.na(curve$se)))
})

test_that("survival curve is monotone non-increasing in time", {
  dt <- sim_constant_hazard(n = 2000L, K = 20L, h = 0.04, seed = 3L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)

  pp_cf <- apply_intervention_pp(fit$pp_data, "A", causatr::static(0))
  haz <- predict_hazard_pp(fit$model, pp_cf)

  curve <- compute_survival_curve(
    pp_data = pp_cf,
    hazards = haz,
    id = "id",
    time = "t",
    times = seq_len(20L),
    intervention_name = "a0"
  )
  expect_true(all(diff(curve$s_hat) <= 0))
})
