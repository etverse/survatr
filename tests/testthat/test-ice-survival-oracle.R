## Track B external point-estimate oracles. Parametric ICE g-computation is
## consistent only when the sequential regressions are correctly specified; the
## survival-tail pseudo-outcome is approximated by a logistic-linear model, so a
## small (~0.01-0.02) finite-model bias vs the truth is expected (lmtp's TMLE
## corrects it; parametric ICE does not). Tolerances reflect that. The tight
## variance pin is the delicatessen oracle; the curve-level pin is the
## forward-simulation truth (test-ice-survival.R).

test_that("Track B matches lmtp_tmle(survival) for static strategies", {
  skip_if_not_installed("lmtp")
  skip_on_cran()

  dat <- sim_ice_feedback(n = 2500L, K = 3L, seed = 21L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ## lmtp's survival fit needs >= 2 outcome nodes; compare at horizons 2..3.
  tt <- 2:3
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = tt,
    type = "survival"
  )
  s1 <- res$estimates[get("intervention") == "a1"][order(get("time"))]$s_hat
  s0 <- res$estimates[get("intervention") == "a0"][order(get("time"))]$s_hat

  truth1 <- true_ice_survival(1, tt, K = 3L)
  truth0 <- true_ice_survival(0, tt, K = 3L)

  lm1 <- oracle_lmtp_survival(dat, a = 1, K = 3L, times = tt)
  lm0 <- oracle_lmtp_survival(dat, a = 0, K = 3L, times = tt)
  skip_if(is.null(lm1) || is.null(lm0), "lmtp oracle did not converge")

  ## lmtp (TMLE, doubly robust) lands near the g-formula truth, validating the
  ## reshape; ICE (parametric) lands near both within its finite-model bias.
  expect_equal(lm1, truth1, tolerance = 0.03)
  expect_equal(lm0, truth0, tolerance = 0.03)
  expect_equal(s1, lm1, tolerance = 0.04)
  expect_equal(s0, lm0, tolerance = 0.04)
})

test_that("Track B matches the forward-sim truth for a dynamic strategy", {
  ## A genuinely longitudinal, non-static strategy: treat iff the current
  ## (time-varying) confounder is positive, A_k = 1{L_k > 0}. Treatment responds
  ## to the evolving covariate, which responds to past treatment -- the
  ## treatment-confounder-feedback case ICE exists for. Forward-sim truth.
  skip_on_cran()
  dat <- sim_ice_feedback(n = 5000L, K = 3L, seed = 27L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  ## causatr calls the rule as `rule(data, orig_trt)`; accept the extra arg.
  rule <- causatr::dynamic(function(data, ...) as.integer(data[["L"]] > 0))
  res <- contrast(
    fit,
    interventions = list(d = rule),
    times = 1:3,
    type = "survival"
  )
  s_d <- res$estimates[order(get("time"))]$s_hat
  truth_d <- true_ice_dynamic(0, 1:3, K = 3L)
  ## Sampling (n = 5000) + parametric-ICE finite-model bias, which accumulates
  ## over horizons (~0.02 at t = 3 on this feedback DGP).
  expect_equal(s_d, truth_d, tolerance = 0.035)
})

test_that("Track B matches gfoRmula::gformula_survival when available", {
  skip_if_not_installed("gfoRmula")
  skip_on_cran()
  ## gfoRmula is the forward-simulation reference named in CHUNK_6. When it is
  ## installed we cross-check the static-strategy risk curve; otherwise this
  ## skips and the forward-sim MC truth (test-ice-survival.R) is the anchor.
  dat <- sim_ice_feedback(n = 2000L, K = 3L, seed = 31L)
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1)),
    times = 1:3,
    type = "risk"
  )
  ice_risk <- res$estimates[order(get("time"))]$risk_hat

  ## gfoRmula reserves `time_name = "t"` and expects a 0-based time column, so
  ## remap t (1..3) -> time0 (0..2). Static interventions take a length-
  ## time_points treatment vector. Wrap defensively so a version-specific API
  ## change skips rather than fails CI.
  gdat <- data.table::copy(data.table::as.data.table(dat))
  gdat[, "time0" := get("t") - 1L][, "t" := NULL]
  gf <- tryCatch(
    gfoRmula::gformula_survival(
      obs_data = gdat,
      id = "id",
      time_name = "time0",
      time_points = 3L,
      covnames = c("L", "A"),
      covtypes = c("normal", "binary"),
      covparams = list(
        covmodels = c(
          L ~ lag1_A + lag1_L,
          A ~ L + lag1_A
        )
      ),
      outcome_name = "Y",
      ymodel = Y ~ A + L,
      intvars = list("A"),
      interventions = list(list(c(gfoRmula::static, rep(1, 3)))),
      int_descript = "Always treat",
      histories = c(gfoRmula::lagged),
      histvars = list(c("A", "L")),
      nsimul = 30000L,
      seed = 1L
    ),
    error = function(e) NULL
  )
  skip_if(is.null(gf), "gfoRmula call did not run with this API")
  ## `g-form risk` at k = 0,1,2 are the intervened risks at horizons 1,2,3.
  res_dt <- data.table::as.data.table(gf$result)
  gf_risk <- res_dt[get("Interv.") == 1][order(get("k"))][["g-form risk"]]
  ## Both are parametric g-computation targeting the same functional via
  ## different algorithms (ICE backward iterated regression vs gfoRmula forward
  ## simulation). On this DGP gfoRmula systematically UNDER-estimates risk by
  ## ~0.02-0.04 (a 3-way check vs the analytic forward-sim truth shows ICE
  ## tracks the truth to ~0.001-0.013 while gfoRmula is the outlier), so the
  ## 0.05 tolerance accommodates gfoRmula's deviation, not ICE's. ICE's tight
  ## point pin is the analytic truth (test-ice-survival.R) and the delicatessen
  ## M-estimator (test-ice-survival-delicatessen.R); this is a ballpark
  ## cross-check against a second g-formula package. Absolute scale: the risks
  ## are small (~0.07-0.29), so `expect_equal`'s relative tolerance would be
  ## inappropriately tight.
  expect_lt(max(abs(ice_risk - as.numeric(gf_risk))), 0.05)
})
