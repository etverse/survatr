## Tests for diagnose() and survatr_diag S3 class.
## DGPs reuse helpers from helper-dgp.R.

## ── Helpers ──────────────────────────────────────────────────────────────────

## Small rectangular person-period fixture for structure / smoke tests.
## No confounding, moderate hazard (no positivity flags expected).
fixture_diag_pp <- function(n = 120L, K = 4L, seed = 20L) {
  set.seed(seed)
  id_col <- rep(seq_len(n), each = K)
  a_col <- rep(stats::rbinom(n, 1L, 0.5), each = K)
  l_col <- rep(stats::rnorm(n), each = K)
  y_col <- stats::rbinom(n * K, 1L, 0.08)
  data.table::data.table(
    id = id_col,
    t = rep(seq_len(K), times = n),
    A = a_col,
    L = l_col,
    Y = y_col
  )
}

## Person-period with a censoring column.
fixture_diag_cens <- function(n = 100L, K = 4L, seed = 21L) {
  pp <- fixture_diag_pp(n, K, seed)
  set.seed(seed + 1L)
  pp[, C := stats::rbinom(.N, 1L, 0.04)]
  pp
}

## Competing-risks person-period fixture.
fixture_diag_cr <- function(n = 100L, K = 4L, seed = 22L) {
  set.seed(seed)
  id_col <- rep(seq_len(n), each = K)
  a_col <- rep(stats::rbinom(n, 1L, 0.5), each = K)
  l_col <- rep(stats::rnorm(n), each = K)
  cause <- integer(n * K)
  for (i in seq_len(n * K)) {
    u <- stats::runif(1L)
    if (u < 0.05) {
      cause[i] <- 1L
    } else if (u < 0.09) {
      cause[i] <- 2L
    }
  }
  data.table::data.table(
    id = id_col,
    t = rep(seq_len(K), times = n),
    A = a_col,
    L = l_col,
    cause = cause
  )
}

## ── Structure ─────────────────────────────────────────────────────────────────

test_that("diagnose.survatr_fit returns a survatr_diag with expected panels", {
  pp <- fixture_diag_pp()
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)

  expect_s3_class(dx, "survatr_diag")
  expect_named(
    dx,
    c("positivity", "balance", "weights", "censoring", "competing"),
    ignore.order = FALSE
  )

  ## Positivity: columns present and row count = number of time periods
  expect_true(data.table::is.data.table(dx$positivity))
  expect_named(
    dx$positivity,
    c("time", "n_at_risk", "h_min", "h_mean", "h_max", "flag_low", "flag_high"),
    ignore.order = FALSE
  )
  expect_equal(nrow(dx$positivity), 4L)

  ## Balance: columns present
  expect_true(data.table::is.data.table(dx$balance))
  expect_true(all(
    c("time", "variable", "smd", "n_a1", "n_a0") %in%
      names(dx$balance)
  ))

  ## weights / censoring / competing NULL for plain gcomp, no censoring column
  expect_null(dx$weights)
  expect_null(dx$censoring)
  expect_null(dx$competing)
})

## ── Positivity panel ─────────────────────────────────────────────────────────

test_that("positivity flags do NOT fire for a moderate-hazard DGP", {
  ## Use sim_constant_hazard (h = 0.08, large n): all model-predicted hazards
  ## should be well within (0.001, 0.999).
  pp <- sim_constant_hazard(n = 300L, K = 4L, h = 0.08, seed = 30L)
  fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)
  expect_false(any(dx$positivity$flag_low))
  expect_false(any(dx$positivity$flag_high))
})

test_that("positivity flag_low fires when predicted hazard < 0.001", {
  ## Very low h → model predicts hazards near 0; glm may not converge due to
  ## near-perfect separation — the convergence warning is expected, not a bug.
  pp <- sim_constant_hazard(n = 400L, K = 3L, h = 1e-4, seed = 31L)
  fit <- withCallingHandlers(
    surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~1),
    warning = function(w) {
      if (grepl("algorithm did not converge", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  dx <- diagnose(fit)
  expect_true(any(dx$positivity$flag_low))
})

## ── Balance panel ─────────────────────────────────────────────────────────────

test_that("balance SMDs are near 0 for a fully randomized DGP", {
  ## `sim_constant_hazard` has no confounders and randomized treatment; the
  ## `L` confounder column is absent, so balance is empty. Use a DGP where L
  ## exists but treatment is independent of L.
  pp <- sim_confounded_survival(n = 400L, K = 4L, seed = 40L, gamma = 0.0)
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)
  expect_true(nrow(dx$balance) > 0L)
  ## With gamma = 0 (no confounding), |SMD| should be small but never exactly 0.
  expect_true(all(abs(dx$balance$smd) < 0.3))
})

test_that("balance SMDs are non-trivial for a confounded DGP", {
  ## With strong confounding (gamma = 1.5) at least some periods have |SMD| > 0.1.
  pp <- sim_confounded_survival(n = 400L, K = 4L, seed = 41L, gamma = 1.5)
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)
  expect_true(any(abs(dx$balance$smd) > 0.1))
})

## ── Weights panel ─────────────────────────────────────────────────────────────

test_that("IPW fit: weights panel is non-NULL and ESS < n_ids", {
  pp <- sim_confounded_survival(n = 200L, K = 4L, seed = 50L, gamma = 0.8)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    estimator = "ipw",
    time_formula = ~1
  )
  dx <- diagnose(fit)
  expect_false(is.null(dx$weights))
  wt <- dx$weights
  expect_equal(wt$n_ids, 200L)
  ## ESS < n for a positively-confounded design (weights are not all 1).
  expect_lt(wt$ess, wt$n_ids)
  ## Sanity: max weight is positive and finite.
  expect_true(is.finite(wt$max_weight) && wt$max_weight > 0)
  ## top5_share in [0, 1]
  expect_true(wt$top5_share >= 0 && wt$top5_share <= 1)
})

test_that("gcomp fit: weights panel is NULL", {
  pp <- fixture_diag_pp()
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)
  expect_null(dx$weights)
})

test_that("ICE fit: weights panel is NULL", {
  ## Use sim_ice_feedback (time-varying treatment) so the once-per-session
  ## survatr_ice_static_treatment inform is NOT consumed here, leaving it
  ## available for the dedicated test in test-ice-survival.R.
  pp <- sim_ice_feedback(n = 60L, K = 3L, seed = 55L)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    estimator = "ice",
    confounders_tv = ~L,
    history = 1L
  )
  dx <- diagnose(fit)
  expect_null(dx$weights)
})

## ── Censoring panel ───────────────────────────────────────────────────────────

test_that("censoring panel is NULL when no censoring column supplied", {
  pp <- fixture_diag_pp()
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  expect_null(diagnose(fit)$censoring)
})

test_that("censoring panel is a data.table with expected columns when C supplied", {
  pp <- fixture_diag_cens()
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "C",
    time_formula = ~1
  )
  dx <- diagnose(fit)
  cens <- dx$censoring
  expect_false(is.null(cens))
  expect_true(data.table::is.data.table(cens))
  expect_true(all(
    c("time", "arm", "n_at_risk", "n_censored", "prop_censored") %in%
      names(cens)
  ))
  ## All proportions in [0, 1]
  expect_true(all(cens$prop_censored >= 0 & cens$prop_censored <= 1))
})

## ── Competing-risks panel ─────────────────────────────────────────────────────

test_that("competing panel: present, correct causes, identity ΣF + S = 1", {
  pp <- fixture_diag_cr()
  fit <- surv_fit(
    pp,
    outcome = "cause",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    time_formula = ~1,
    competing = "cause"
  )
  dx <- diagnose(fit)

  expect_false(is.null(dx$competing))
  cr <- dx$competing
  expect_true(data.table::is.data.table(cr))
  expect_true(all(c("cause", "time", "cif_hat") %in% names(cr)))
  ## Both causes present in the panel.
  expect_true(all(c(1L, 2L) %in% cr$cause))

  ## Identity check: ΣF^(j)(t) + S(t) = 1 to floating-point precision.
  id_chk <- attr(cr, "identity_check")
  expect_false(is.null(id_chk))
  expect_lt(max(id_chk$abs_dev, na.rm = TRUE), 1e-6)

  ## Caveat attribute is a non-empty string.
  cav <- attr(cr, "caveat")
  expect_true(is.character(cav) && nchar(cav) > 0L)
})

test_that("competing panel is NULL for a single-event fit", {
  pp <- fixture_diag_pp()
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  expect_null(diagnose(fit)$competing)
})

## ── print.survatr_diag ────────────────────────────────────────────────────────

test_that("print.survatr_diag runs without error and returns invisibly", {
  pp <- fixture_diag_pp()
  fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
  dx <- diagnose(fit)
  out <- expect_no_error(capture.output(print(dx)))
  expect_true(any(grepl("survatr_diag", out)))
})

test_that("print.survatr_diag shows weights line for IPW", {
  pp <- sim_confounded_survival(n = 100L, K = 3L, seed = 60L, gamma = 0.8)
  fit <- surv_fit(
    pp,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    estimator = "ipw",
    time_formula = ~1
  )
  out <- capture.output(print(diagnose(fit)))
  expect_true(any(grepl("Weights:", out)))
})

test_that("print.survatr_diag shows competing line for CR fit", {
  pp <- fixture_diag_cr()
  fit <- surv_fit(
    pp,
    outcome = "cause",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    time_formula = ~1,
    competing = "cause"
  )
  out <- capture.output(print(diagnose(fit)))
  expect_true(any(grepl("Competing:", out)))
})
