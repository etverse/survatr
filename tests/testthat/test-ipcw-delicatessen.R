## IPCW three-block stacked-EE sandwich vs an INDEPENDENT analytic oracle:
## delicatessen (Zivich 2022) stacked M-estimation in Python. Both
## implementations operate on the same committed person-period fixture with
## informative censoring (delta_cens = 0.8), so point-estimate agreement is to
## solver tolerance and SE agreement is expected within ~2% (independent code
## paths, same three-block EE structure).
##
## This test pins:
##   (a) the two-block treatment-model bread projection (Block 1)
##   (b) the IPCW censoring-model bread projection (Block 2)
##   (c) the combined three-block running-product delta chain
##
## Reference values generated offline by
## `data-raw/delicatessen_ipcw_survival.py` (run with causatr's delicatessen
## venv) and committed as `fixtures/python/ipcw_survival_delicatessen.csv`. See
## `fixtures/python/README.md` to regenerate.

test_that("IPCW three-block sandwich matches the delicatessen M-estimation oracle", {
  data_fx <- test_path("fixtures", "python", "ipcw_survival_data.csv")
  ref_fx <- test_path("fixtures", "python", "ipcw_survival_delicatessen.csv")
  skip_if(!file.exists(ref_fx), "delicatessen IPCW reference fixture not generated")

  dt <- data.table::fread(data_fx)
  ref <- data.table::fread(ref_fx)
  tt <- sort(unique(ref$time))

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
  sv <- contrast(
    fit,
    interventions = ivs,
    times = tt,
    type = "survival",
    ci_method = "sandwich"
  )$estimates
  rd <- contrast(
    fit,
    interventions = ivs,
    times = tt,
    type = "risk_difference",
    reference = "a1",
    ci_method = "sandwich"
  )$contrasts

  ref_s1 <- ref[get("quantity") == "S1"][order(get("time"))]
  ref_s0 <- ref[get("quantity") == "S0"][order(get("time"))]
  ref_rd <- ref[get("quantity") == "RD"][order(get("time"))]
  s1 <- sv[get("intervention") == "a1"][order(get("time"))]
  s0 <- sv[get("intervention") == "a0"][order(get("time"))]
  rd <- rd[order(get("time"))]

  ## Point estimates: identical data + MSM design, agree to solver tolerance.
  expect_equal(s1$s_hat, ref_s1$estimate, tolerance = 1e-4)
  expect_equal(s0$s_hat, ref_s0$estimate, tolerance = 1e-4)
  expect_equal(rd$estimate, ref_rd$estimate, tolerance = 1e-4)

  ## SEs: survatr's analytic three-block cross-time delta chain vs delicatessen's
  ## numerical stacked-EE sandwich. The censoring-model numerator is omitted from
  ## the IF in both implementations (conservative, consistent with Robins 1999 for
  ## stabilized weights). Agreement within 2% pins the cross-derivative
  ## A_beta_gamma scaling and the n_ids/n_cens_fit bread correction.
  expect_equal(s1$se, ref_s1$se, tolerance = 0.02)
  expect_equal(s0$se, ref_s0$se, tolerance = 0.02)
  expect_equal(rd$se, ref_rd$se, tolerance = 0.02)
})
