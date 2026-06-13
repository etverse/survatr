## RMTL (restricted mean time lost) is a smooth functional of the survival
## curve already produced by chunks 2/3: RMTL(t) = t - RMST(t). The tests pin
## the two defining identities -- the point identity `RMTL = t - RMST` and the
## variance identity `Var(RMTL) = Var(RMST)` (the constant t drops out of the
## delta gradient) -- so RMTL inherits the RMST oracle pins (survRM2,
## delicatessen) transitively rather than re-deriving an independent oracle.

test_that("RMTL = t - RMST identity holds pointwise (sandwich path)", {
  dt <- sim_constant_hazard(n = 400L, K = 8L, h = 0.06, seed = 401L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  times <- 1:8

  r_rmst <- contrast(fit, ivs, times, type = "rmst", ci_method = "sandwich")
  r_rmtl <- contrast(fit, ivs, times, type = "rmtl", ci_method = "sandwich")

  for (iv in c("a1", "a0")) {
    e_rmst <- r_rmst$estimates[get("intervention") == iv][order(time)]
    e_rmtl <- r_rmtl$estimates[get("intervention") == iv][order(time)]
    ## Exact algebraic identity: machine precision, not simulation tolerance.
    expect_equal(e_rmtl$rmtl_hat, times - e_rmst$rmst_hat, tolerance = 1e-12)
  }
})

test_that("Var(RMTL) = Var(RMST): the SE is identical", {
  dt <- sim_constant_hazard(n = 400L, K = 8L, h = 0.06, seed = 409L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  times <- 1:8

  se_rmst <- contrast(
    fit,
    ivs,
    times,
    type = "rmst",
    ci_method = "sandwich"
  )$estimates[get("intervention") == "a1"][order(time)]$se
  se_rmtl <- contrast(
    fit,
    ivs,
    times,
    type = "rmtl",
    ci_method = "sandwich"
  )$estimates[get("intervention") == "a1"][order(time)]$se
  ## The constant restriction time has zero gradient, so the trapezoidal
  ## quadratic form -- and hence the SE -- is the same as RMST's.
  expect_equal(se_rmtl, se_rmst, tolerance = 1e-12)
  expect_true(all(se_rmtl >= 0))
})

test_that("rmtl_difference is the negative of rmst_difference (point + SE)", {
  dt <- sim_confounded_survival(n = 600L, K = 6L, seed = 415L)
  fit <- surv_fit(dt, "Y", "A", ~L, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  times <- 1:6

  d_rmst <- contrast(
    fit,
    ivs,
    times,
    type = "rmst_difference",
    reference = "a0",
    ci_method = "sandwich"
  )$contrasts[order(time)]
  d_rmtl <- contrast(
    fit,
    ivs,
    times,
    type = "rmtl_difference",
    reference = "a0",
    ci_method = "sandwich"
  )$contrasts[order(time)]

  ## RMTL(t) = t - RMST(t) => RMTL_a1 - RMTL_a0 = -(RMST_a1 - RMST_a0).
  expect_equal(d_rmtl$estimate, -d_rmst$estimate, tolerance = 1e-12)
  ## Variance ignores the sign, so the SE and CI half-width match exactly.
  expect_equal(d_rmtl$se, d_rmst$se, tolerance = 1e-12)
})

test_that("sandwich and bootstrap SE agree for RMTL", {
  dt <- sim_confounded_survival(n = 800L, K = 6L, seed = 421L)
  fit <- surv_fit(dt, "Y", "A", ~L, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  times <- 1:6

  se_sand <- contrast(
    fit,
    ivs,
    times,
    type = "rmtl",
    ci_method = "sandwich"
  )$estimates[get("intervention") == "a1"][order(time)]$se
  se_boot <- contrast(
    fit,
    ivs,
    times,
    type = "rmtl",
    ci_method = "bootstrap",
    n_boot = 150L,
    seed = 99L
  )$estimates[get("intervention") == "a1"][order(time)]$se

  ## Two SE estimators of the same quantity: ratio within (0.7, 1.4) at the
  ## last (most variable) time point.
  ratio <- se_boot[length(times)] / se_sand[length(times)]
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})

test_that("RMTL is wired across estimators (IPW, IPCW, ICE inheritance)", {
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  ## IPW: stabilized weighted MSM -> same survival IF -> RMTL for free.
  dt_ipw <- sim_confounded_survival(n = 500L, K = 5L, seed = 431L)
  fit_ipw <- surv_fit(
    dt_ipw,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  r_ipw <- contrast(fit_ipw, ivs, 1:5, type = "rmtl", ci_method = "sandwich")
  expect_true(all(is.finite(r_ipw$estimates$rmtl_hat)))
  expect_true(all(r_ipw$estimates$se >= 0))

  ## IPCW: per-period censoring weights -> third sandwich block -> free RMTL.
  dt_ipcw <- sim_informative_censoring(n = 400L, K = 5L, seed = 433L)
  fit_ipcw <- surv_fit(
    dt_ipcw,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "C",
    time_formula = ~ factor(t),
    estimator = "ipw",
    ipcw = ~L
  )
  r_ipcw <- contrast(fit_ipcw, ivs, 1:5, type = "rmtl", ci_method = "sandwich")
  expect_true(all(is.finite(r_ipcw$estimates$rmtl_hat)))
  expect_true(all(r_ipcw$estimates$se >= 0))

  ## Track B (ICE): the shared fill_sandwich_ses() consumes the ICE IF matrix,
  ## so RMTL is available with no Track-B-specific code.
  dt_ice <- sim_ice_feedback(n = 500L, K = 4L, seed = 437L)
  fit_ice <- surv_fit(
    dt_ice,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  r_ice <- contrast(fit_ice, ivs, 1:4, type = "rmtl", ci_method = "sandwich")
  expect_true(all(is.finite(r_ice$estimates$rmtl_hat)))
  expect_true(all(r_ice$estimates$se >= 0))
})

test_that("RMTL flows through tidy() and plot()", {
  dt <- sim_constant_hazard(n = 600L, K = 5L, h = 0.07, seed = 441L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  res <- contrast(
    fit,
    ivs,
    1:5,
    type = "rmtl_difference",
    ci_method = "sandwich"
  )
  td <- tidy(res)
  expect_true("rmtl_hat" %in% td$estimand)
  expect_true("rmtl_difference" %in% td$estimand)

  pf <- tempfile(fileext = ".png")
  grDevices::png(pf)
  plot(res)
  grDevices::dev.off()
  expect_true(file.exists(pf))
})
