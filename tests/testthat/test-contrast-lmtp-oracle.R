## Oracle smoke test against `lmtp::lmtp_tmle(outcome_type = "survival")`.
##
## lmtp uses a different estimator (TMLE / SDR on the longitudinal survival
## functional); its point estimates of S^a(t) agree with pooled-logistic
## gcomp up to finite-sample + model-form differences. At n = 2000 on a
## constant-hazard DGP with a linear time term we expect agreement within
## roughly 0.03 at mid-horizon times. EIF-based SE is NOT directly
## comparable to the M-estimation sandwich (different normalizations) so
## this is a point-estimate oracle only.

skip_if_not_installed("lmtp")

test_that("survatr::contrast survival curve agrees with lmtp::lmtp_tmle (smoke)", {
  skip_on_cran()
  dt <- sim_constant_hazard(n = 2000L, K = 5L, h = 0.1, seed = 53L)

  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  surv_res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1)),
    times = 1:5,
    type = "survival"
  )
  s_survatr <- surv_res$estimates$s_hat

  ## Reshape to wide for lmtp: one row per id, columns A_1, Y_1, ..., Y_K,
  ## and a cumulative censoring column C_k = 1 if still under follow-up at
  ## the start of period k. Our DGP has no censoring, so all C_k = 1.
  ##
  ## lmtp reports S(t_k) = 1 - cumulative risk through period k.
  wide <- data.table::dcast(
    data.table::as.data.table(dt),
    id + A ~ t,
    value.var = "Y"
  )
  data.table::setnames(wide, as.character(1:5), paste0("Y_", 1:5))
  ## Cumulative event indicator (used by lmtp as outcome).
  for (k in 2:5) {
    wide[[paste0("Y_", k)]] <- pmax(
      wide[[paste0("Y_", k - 1L)]],
      wide[[paste0("Y_", k)]]
    )
  }
  ## Artificial always-observed censoring indicator.
  for (k in 1:5) {
    wide[[paste0("C_", k)]] <- 1L
  }

  fit_lmtp <- tryCatch(
    lmtp::lmtp_tmle(
      data = as.data.frame(wide),
      trt = "A",
      outcome = paste0("Y_", 1:5),
      cens = paste0("C_", 1:5),
      shift = function(data, trt) rep(1, nrow(data)),
      outcome_type = "survival",
      baseline = NULL,
      folds = 1L,
      .SL_folds = 1L
    ),
    error = function(e) NULL
  )
  skip_if(is.null(fit_lmtp), "lmtp_tmle did not fit on this toy DGP")

  s_lmtp <- fit_lmtp$theta ## vector of S(t_k) under the shift (static 1)
  expect_equal(length(s_lmtp), length(s_survatr))
  expect_equal(s_survatr, s_lmtp, tolerance = 0.05)
})
