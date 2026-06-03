## Track A IPW (weighted hazard MSM) -- point estimate, structural invariants,
## and the IPW-specific rejection surface. Variance is in
## `test-ipw-survival-sandwich.R` / `test-ipw-survival-bootstrap.R`.

test_that("IPW survival curve matches lmtp::lmtp_tmle on a confounded DGP", {
  skip_on_cran()
  skip_if_not_installed("lmtp")
  ## Confounded point treatment: A depends on L which also drives the hazard.
  ## lmtp's TMLE survival functional is the external point-estimate oracle.
  dt <- sim_confounded_survival(n = 2000L, K = 5L, seed = 3L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "survival"
  )
  s1 <- res$estimates[get("intervention") == "a1"]$s_hat
  s0 <- res$estimates[get("intervention") == "a0"]$s_hat

  ## Compare at mid-horizon and the endpoint (two lmtp fits per arm keeps the
  ## test affordable; the endpoint is the most accumulated point of the curve).
  tt <- c(3L, 5L)
  o1 <- lmtp_ipw_survival_oracle(dt, 1, times = tt)
  o0 <- lmtp_ipw_survival_oracle(dt, 0, times = tt)
  skip_if(is.null(o1) || is.null(o0), "lmtp_tmle did not fit the toy DGP")

  expect_equal(s1[tt], o1, tolerance = 0.05)
  expect_equal(s0[tt], o0, tolerance = 0.05)
})

test_that("single-period IPW reduces to causatr's scalar point IPW", {
  skip_on_cran()
  ## Degenerate-to-scalar check: with one period the survival MSM is saturated
  ## in binary A, so survatr's stabilized-weight IPW risk equals causatr's
  ## Hajek point-IPW E[Y^a]. This pins the engine boundary.
  set.seed(5)
  n <- 3000L
  L <- rnorm(n)
  A <- rbinom(n, 1L, plogis(0.8 * L))
  Y <- rbinom(n, 1L, plogis(-0.5 - 0.7 * A + 0.6 * L))
  pp <- data.table::data.table(id = seq_len(n), t = 1L, A = A, L = L, Y = Y)

  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~1,
    estimator = "ipw"
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1,
    type = "risk"
  )
  r1 <- res$estimates[get("intervention") == "a1"]$risk_hat
  r0 <- res$estimates[get("intervention") == "a0"]$risk_hat

  scal <- data.frame(A = A, L = L, Y = Y)
  cf <- causatr::causat(
    scal,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    estimator = "ipw",
    model_fn = stats::glm,
    propensity_model_fn = stats::glm
  )
  ct <- causatr::contrast(
    cf,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0))
  )
  mu1 <- ct$estimates[get("intervention") == "a1"]$estimate
  mu0 <- ct$estimates[get("intervention") == "a0"]$estimate

  expect_equal(r1, mu1, tolerance = 1e-2)
  expect_equal(r0, mu0, tolerance = 1e-2)
})

test_that("unconfounded IPW (weights ~ 1) matches the unadjusted gcomp curve", {
  ## With A independent of L the stabilized weights collapse to ~1, so the
  ## weighted marginal MSM coincides with the unadjusted gcomp curve.
  dt <- sim_confounded_survival(n = 3000L, K = 5L, seed = 8L, gamma = 0)
  fit_ipw <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  ## Weights are near 1 and tight.
  expect_equal(mean(fit_ipw$weights), 1, tolerance = 0.05)
  expect_lt(stats::sd(fit_ipw$weights), 0.25)

  fit_gc <- surv_fit(
    dt,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "gcomp"
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  s_ipw <- contrast(
    fit_ipw,
    interventions = ivs,
    times = 1:5,
    type = "survival"
  )$estimates$s_hat
  s_gc <- contrast(
    fit_gc,
    interventions = ivs,
    times = 1:5,
    type = "survival"
  )$estimates$s_hat
  expect_equal(s_ipw, s_gc, tolerance = 0.02)
})

test_that("the stabilized weight is constant within id (broadcast invariant)", {
  dt <- sim_confounded_survival(n = 200L, K = 4L, seed = 2L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  w_dt <- data.table::data.table(id = dt$id, w = fit$weights)
  n_uniq <- w_dt[, list(nu = data.table::uniqueN(w)), by = "id"]$nu
  expect_true(all(n_uniq == 1L))
  ## fit$weights aligns with the person-period rows.
  expect_length(fit$weights, nrow(dt))
})

test_that("trim winsorizes the per-id weights at a fixed cutoff", {
  dt <- sim_confounded_survival(n = 1000L, K = 4L, seed = 4L, gamma = 1.2)
  fit_raw <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  fit_trim <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw",
    trim = 0.9
  )
  expect_true(is.finite(fit_trim$trim_threshold))
  ## No weight exceeds the resolved cutoff, and trimming lowers the maximum.
  expect_lte(max(fit_trim$weights), fit_trim$trim_threshold + 1e-8)
  expect_lt(max(fit_trim$weights), max(fit_raw$weights))
  ## Untrimmed fit records NA threshold.
  expect_true(is.na(fit_raw$trim_threshold))
})

test_that("IPW is robust to a collinear confounder (aliased coef dropped)", {
  dt <- sim_confounded_survival(n = 800L, K = 4L, seed = 6L, gamma = 0.8)
  dt[, L2 := 2 * L] # perfectly collinear with L
  ## causatr drops the aliased propensity coefficient with a warning; the
  ## fit and curve must still succeed and be finite.
  fit <- suppressWarnings(surv_fit(
    dt,
    "Y",
    "A",
    ~ L + L2,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  ))
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:4,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(is.finite(res$contrasts$se)))
})

test_that("IPW rejects unsupported / invalid combinations (classed)", {
  base <- sim_confounded_survival(n = 300L, K = 4L, seed = 1L, gamma = 0.8)

  ## Continuous treatment -> deferred extended-types path.
  dt_cont <- data.table::copy(base)
  set.seed(1)
  a_cont <- rnorm(length(unique(dt_cont$id)))
  dt_cont[, A := a_cont[match(id, sort(unique(id)))]]
  expect_error(
    surv_fit(dt_cont, "Y", "A", ~L, "id", "t", estimator = "ipw"),
    class = "survatr_ipw_treatment_unsupported"
  )

  ## Time-varying treatment within id.
  dt_tv <- data.table::copy(base)
  dt_tv[, A := as.integer(t %% 2L)] # varies within id
  expect_error(
    surv_fit(dt_tv, "Y", "A", ~L, "id", "t", estimator = "ipw"),
    class = "survatr_ipw_time_varying_treatment"
  )

  ## External weights composed with IPW.
  expect_error(
    surv_fit(
      base,
      "Y",
      "A",
      ~L,
      "id",
      "t",
      estimator = "ipw",
      weights = rep(1, nrow(base))
    ),
    class = "survatr_ipw_external_weights"
  )
})

test_that("ipsi() on an IPW fit is rejected as deferred", {
  dt <- sim_confounded_survival(n = 300L, K = 4L, seed = 1L, gamma = 0.8)
  fit <- surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    time_formula = ~ factor(t),
    estimator = "ipw"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = causatr::static(1), shift = causatr::ipsi(2)),
      times = 1:4,
      type = "survival"
    ),
    class = "survatr_ipw_ipsi_deferred"
  )
})

test_that("IPW rejects degenerate treatment and NA confounders (classed)", {
  base <- sim_confounded_survival(n = 300L, K = 4L, seed = 1L, gamma = 0.8)

  ## Constant treatment across all ids: no treated/untreated contrast, so the
  ## inverse-probability weights are undefined (positivity violation). Must give
  ## a clear no-variation message, not causatr's "gaussian/continuous" one.
  d_const <- data.table::copy(base)
  d_const[, A := 1L]
  expect_error(
    surv_fit(d_const, "Y", "A", ~L, "id", "t", estimator = "ipw"),
    class = "survatr_ipw_no_treatment_variation"
  )

  ## NA in a confounder is rejected upfront. survatr does not delegate missing-
  ## data handling to causatr's row-dropping: the IPW treatment model's
  ## predictors flow through `check_no_na_in_predictors()` like the gcomp path,
  ## so the influence-function row alignment is preserved.
  d_na <- data.table::copy(base)
  d_na[1L, L := NA_real_]
  expect_error(
    surv_fit(d_na, "Y", "A", ~L, "id", "t", estimator = "ipw"),
    class = "survatr_na_in_predictors"
  )
})
