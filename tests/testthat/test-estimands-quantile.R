## The survival quantile is a smooth functional of the curve chunks 2/3
## already compute: tau_q = inf{t : S(t) <= 1 - q}. On a constant-hazard DGP
## the discrete survival is S(t) = (1 - h)^t, so the q-quantile has the closed
## form tau_q = -log(1 - q) / (-log(1 - h)) (median = log(2) / (-log(1 - h))).
## The estimator returns the linear-interpolation crossing of the FITTED curve,
## which converges to that closed form as h_hat -> h. n is large here so the
## fitted hazard is close to nominal; the residual is a ~0.3% downward
## discretization bias (linear interpolation of a convex curve, unit spacing).

test_that("median and quantiles match the constant-hazard closed form", {
  h <- 0.05
  ## Large n: at n = 4000 the per-period MLE h_hat can sit ~5% off nominal,
  ## shifting the median by the same factor (verified: estimator returns the
  ## exact linear-interp crossing of the fitted curve at every n). n = 20000
  ## pins h_hat to within ~0.3% so the comparison to the NOMINAL rate is tight.
  dt <- sim_constant_hazard(n = 20000L, K = 60L, h = h, seed = 501L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  ivs <- list(a0 = causatr::static(0))
  lambda <- -log(1 - h) ## discrete-time rate per period

  ## Median (q = 0.5).
  med <- contrast(fit, ivs, 1:60, type = "quantile", q = 0.5)$estimates$tau_hat
  expect_equal(med, log(2) / lambda, tolerance = 0.02)

  ## A vector of quantiles, each against -log(1 - q) / lambda.
  qs <- c(0.25, 0.5, 0.75)
  res <- contrast(fit, ivs, 1:60, type = "quantile", q = qs)
  truth <- -log(1 - qs) / lambda
  expect_equal(res$estimates[order(q)]$tau_hat, truth, tolerance = 0.03)
})

test_that("quantile sandwich and bootstrap SE agree", {
  dt <- sim_constant_hazard(n = 3000L, K = 40L, h = 0.06, seed = 511L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  se_sand <- contrast(
    fit,
    ivs,
    1:40,
    type = "quantile",
    q = 0.5,
    ci_method = "sandwich"
  )$estimates[get("intervention") == "a1"]$se
  se_boot <- contrast(
    fit,
    ivs,
    1:40,
    type = "quantile",
    q = 0.5,
    ci_method = "bootstrap",
    n_boot = 200L,
    seed = 7L
  )$estimates[get("intervention") == "a1"]$se

  ratio <- se_boot / se_sand
  expect_gt(ratio, 0.7)
  expect_lt(ratio, 1.4)
})

test_that("survatr median agrees with survival::survfit on the unadjusted KM", {
  skip_if_not_installed("survival")
  ## Pooled-logistic median ~ KM median for small per-period hazard. Build the
  ## one-row-per-id (time, status) table the KM needs from the PP grid.
  dt <- sim_constant_hazard(n = 3000L, K = 50L, h = 0.05, seed = 521L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  med_survatr <- contrast(
    fit,
    list(a0 = causatr::static(0)),
    1:50,
    type = "quantile",
    q = 0.5
  )$estimates$tau_hat

  ## Per-id event time = first event period (or censored at K).
  per_id <- dt[,
    list(
      time = if (any(Y == 1L)) min(t[Y == 1L]) else max(t),
      status = as.integer(any(Y == 1L))
    ),
    by = "id"
  ]
  km <- survival::survfit(survival::Surv(time, status) ~ 1, data = per_id)
  med_km <- summary(km)$table[["median"]]
  ## KM median is a step-function jump; pooled-logistic interpolates. Agree to
  ## within a couple of periods.
  expect_equal(med_survatr, med_km, tolerance = 0.1)
})

test_that("median difference contrast is well-formed", {
  dt <- sim_confounded_survival(n = 2000L, K = 12L, seed = 531L)
  fit <- surv_fit(dt, "Y", "A", ~L, "id", "t", time_formula = ~ factor(t))
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  res <- contrast(
    fit,
    ivs,
    1:12,
    type = "quantile",
    q = 0.5,
    reference = "a0",
    ci_method = "sandwich"
  )
  ## One contrast row (a1 vs a0) at q = 0.5; estimate = tau(a1) - tau(a0).
  expect_equal(nrow(res$contrasts), 1L)
  tau_a1 <- res$estimates[get("intervention") == "a1"]$tau_hat
  tau_a0 <- res$estimates[get("intervention") == "a0"]$tau_hat
  expect_equal(res$contrasts$estimate, tau_a1 - tau_a0, tolerance = 1e-10)
  expect_true(is.finite(res$contrasts$se))
})

test_that("a single intervention is accepted (no pairwise contrast required)", {
  dt <- sim_constant_hazard(n = 800L, K = 30L, h = 0.06, seed = 541L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  ## Unlike the strict *_difference types, a lone median is a valid request.
  res <- contrast(
    fit,
    list(a0 = causatr::static(0)),
    1:30,
    type = "quantile",
    q = 0.5,
    ci_method = "sandwich"
  )
  expect_equal(nrow(res$estimates), 1L)
  expect_equal(nrow(res$contrasts), 0L)
  expect_true(is.finite(res$estimates$se))
})

test_that("single-intervention quantile bootstrap does not crash", {
  ## The quantile is the only contrast-kind estimand that bypasses the
  ## >=2-intervention guard, so a single arm reaches the bootstrap layout
  ## code. `setdiff(iv, ref)` is then empty and `paste0()` used to recycle it
  ## into a phantom contrast column, crashing every replicate ("replacement has
  ## length zero"). The bootstrap must run and return a finite SE ~ sandwich.
  dt <- sim_constant_hazard(n = 1500L, K = 30L, h = 0.06, seed = 571L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  arm <- list(a0 = causatr::static(0))
  se_sand <- contrast(
    fit,
    arm,
    1:30,
    type = "quantile",
    q = 0.5,
    ci_method = "sandwich"
  )$estimates$se
  res_boot <- contrast(
    fit,
    arm,
    1:30,
    type = "quantile",
    q = 0.5,
    ci_method = "bootstrap",
    n_boot = 150L,
    seed = 2L
  )
  expect_equal(nrow(res_boot$estimates), 1L)
  expect_equal(nrow(res_boot$contrasts), 0L)
  expect_true(is.finite(res_boot$estimates$se))
  ratio <- res_boot$estimates$se / se_sand
  expect_gt(ratio, 0.6)
  expect_lt(ratio, 1.5)
})

test_that("q-indexed quantile flows through plot() but not forrest()", {
  dt <- sim_constant_hazard(n = 600L, K = 30L, h = 0.06, seed = 573L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  ## plot() default (auto) on a single-arm quantile falls back to the curve
  ## view rather than aborting on the empty contrasts table.
  solo <- contrast(
    fit,
    list(a0 = causatr::static(0)),
    1:30,
    type = "quantile",
    q = c(0.25, 0.5, 0.75),
    ci_method = "sandwich"
  )
  pf <- tempfile(fileext = ".png")
  grDevices::png(pf)
  plot(solo)
  grDevices::dev.off()
  expect_true(file.exists(pf))

  ## forrest() is a time-slice and is meaningless for the q-indexed quantile;
  ## it aborts with a classed error rather than a raw data.table `==` error.
  pair <- contrast(
    fit,
    ivs,
    1:30,
    type = "quantile",
    q = 0.5,
    reference = "a0",
    ci_method = "sandwich"
  )
  expect_error(forrest(pair, t_ref = 30), class = "survatr_forrest_wrong_type")
})

test_that("quantile is wired across estimators (IPW / ICE / competing risks)", {
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))

  ## IPW: stabilized weighted MSM -> same survival IF -> quantile for free.
  dt_ipw <- sim_confounded_survival(n = 800L, K = 12L, seed = 551L)
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
  r_ipw <- contrast(
    fit_ipw,
    ivs,
    1:12,
    type = "quantile",
    q = 0.5,
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(r_ipw$estimates$tau_hat)))
  expect_true(all(r_ipw$estimates$se >= 0))

  ## Track B (ICE): the shared quantile assembly consumes the ICE IF matrix.
  ## q = 0.3 is reachable on this high-survival 4-period feedback DGP.
  dt_ice <- sim_ice_feedback(n = 1500L, K = 4L, seed = 553L)
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
  r_ice <- contrast(
    fit_ice,
    ivs,
    1:4,
    type = "quantile",
    q = 0.3,
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(r_ice$estimates$tau_hat)))
  expect_true(all(r_ice$estimates$se >= 0))

  ## Competing risks: the all-cause survival quantile (no cause dimension).
  dt_cr <- sim_two_cause_constant_hazard(n = 1500L, K = 25L, seed = 557L)
  fit_cr <- surv_fit(
    dt_cr,
    "event",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    competing = "event"
  )
  r_cr <- contrast(
    fit_cr,
    ivs,
    1:25,
    type = "quantile",
    q = 0.5,
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(r_cr$estimates$tau_hat)))
  expect_false("cause" %in% names(r_cr$estimates)) ## all-cause: no cause dim
})

test_that("unreached quantile and bad q abort with distinct classes", {
  dt <- sim_constant_hazard(n = 600L, K = 6L, h = 0.05, seed = 561L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  ivs <- list(a0 = causatr::static(0))
  ## The curve does not drop to 1 - 0.95 = 0.05 within 6 periods.
  expect_error(
    contrast(fit, ivs, 1:6, type = "quantile", q = 0.95),
    class = "survatr_quantile_unreached"
  )
  ## q must be strictly inside (0, 1).
  expect_error(
    contrast(fit, ivs, 1:6, type = "quantile", q = 1.5),
    class = "survatr_bad_q"
  )
  expect_error(
    contrast(fit, ivs, 1:6, type = "quantile", q = 0),
    class = "survatr_bad_q"
  )
})
