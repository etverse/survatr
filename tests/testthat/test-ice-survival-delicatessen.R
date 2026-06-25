## longitudinal ICE-hazard sandwich vs an INDEPENDENT analytic oracle:
## delicatessen (Zivich 2022) stacked M-estimation in Python. Both
## implementations run the SAME backward survival-tail ICE estimator on the
## same committed person-period fixture, and delicatessen forms the sandwich
## numerically from the stacked estimating equations. Agreement to solver
## tolerance is a far tighter check than the bootstrap and, crucially, an
## independent confirmation of the (1 - D_k) failure carry-forward factor in
## survatr's survival influence-function chain (which causatr's plain ICE chain
## does not have).
##
## Reference values are generated offline by
## `data-raw/delicatessen_ice_survival.py` (run with causatr's delicatessen
## venv) and committed as `fixtures/python/ice_survival_delicatessen.csv`. See
## `fixtures/python/README.md` to regenerate.

test_that("longitudinal ICE-hazard ICE sandwich matches the delicatessen M-estimation oracle", {
  data_fx <- test_path("fixtures", "python", "ice_survival_data.csv")
  ref_fx <- test_path("fixtures", "python", "ice_survival_delicatessen.csv")
  skip_if(!file.exists(ref_fx), "delicatessen reference fixture not generated")

  dt <- data.table::fread(data_fx)
  ref <- data.table::fread(ref_fx)
  tt <- sort(unique(ref$time))

  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  sv <- contrast(
    fit,
    interventions = ivs,
    times = tt,
    type = "survival",
    ci_method = "sandwich"
  )$estimates
  ## reference = a1 makes the "a0 vs a1" contrast risk^0 - risk^1 = S^1 - S^0,
  ## matching how the Python script defines RD = mu0 - mu1.
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

  ## Point estimates: same data + same ICE recursion, two different solvers ->
  ## agree to solver tolerance.
  expect_equal(s1$s_hat, ref_s1$estimate, tolerance = 1e-4)
  expect_equal(s0$s_hat, ref_s0$estimate, tolerance = 1e-4)
  expect_equal(rd$estimate, ref_rd$estimate, tolerance = 1e-4)

  ## Sandwich SEs: survatr's analytic survival IF chain (Channel 1 + the
  ## (1 - D_k)-corrected nuisance cascade) vs delicatessen's numerical
  ## stacked-EE sandwich. Independent code paths; they agree to ~1e-4, so 1%
  ## is a comfortably tight pin that still guards the failure-carry-forward
  ## factor (dropping it inflates these SEs by tens of percent at later t).
  expect_equal(s1$se, ref_s1$se, tolerance = 0.01)
  expect_equal(s0$se, ref_s0$se, tolerance = 0.01)
  expect_equal(rd$se, ref_rd$se, tolerance = 0.01)
})
