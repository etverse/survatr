## `fixture_small_pp()` has only 10 at-risk rows with 2 events, so pooled
## logistic with a continuous confounder and a linear time term is near
## separation — `glm.fit()` emits convergence / fitted-probability warnings
## that are expected on a fixture this tiny. The structural assertions below
## are the point of these tests; the truth-based family-switch oracle lives
## in test-surv_fit-family-oracle.R on a large DGP.

test_that("surv_fit returns a survatr_fit with expected slots", {
  dt <- fixture_small_pp()
  fit <- suppressWarnings(surv_fit(
    data = dt,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    censoring = "cens",
    time_formula = ~t
  ))
  expect_s3_class(fit, "survatr_fit")
  expect_s3_class(fit$model, "glm")
  expect_equal(fit$track, "A")
  expect_equal(fit$estimator, "gcomp")
  expect_equal(fit$family, "binomial")
  expect_equal(fit$outcome, "Y")
  expect_equal(fit$treatment, "A")
  expect_equal(fit$id, "id")
  expect_equal(fit$time, "t")
  expect_equal(fit$censoring, "cens")
  expect_equal(fit$time_grid, 1:3)
  expect_equal(fit$n_total, 15L)
  expect_equal(fit$n_fit, 10L)
})

test_that("surv_fit strips internal bookkeeping columns from pp_data", {
  dt <- fixture_small_pp()
  fit <- suppressWarnings(surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "cens",
    time_formula = ~t
  ))
  expect_false(".survatr_prev_event" %in% names(fit$pp_data))
  expect_false(".survatr_prev_cens" %in% names(fit$pp_data))
})

test_that("surv_fit does not mutate the caller's data", {
  dt <- fixture_small_pp()
  snap_names <- names(dt)
  snap_nrow <- nrow(dt)
  suppressWarnings(surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "cens",
    time_formula = ~t
  ))
  expect_identical(names(dt), snap_names)
  expect_equal(nrow(dt), snap_nrow)
})

test_that("surv_fit rejects matching estimator with pointer error", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "matching"),
    class = "survatr_matching_rejected"
  )
  expect_snapshot(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "matching"),
    error = TRUE
  )
})

test_that("surv_fit rejects ipw / ice / unknown estimators", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "ipw"),
    class = "survatr_bad_estimator"
  )
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "ice"),
    class = "survatr_bad_estimator"
  )
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "bogus"),
    class = "survatr_bad_estimator"
  )
})

test_that("surv_fit rejects non-NULL competing argument", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(
      dt,
      "Y",
      "A",
      ~L,
      "id",
      "t",
      competing = "event_type"
    ),
    class = "survatr_competing_misuse"
  )
})

test_that("surv_fit rejects na.exclude via ... gate", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(
      dt,
      "Y",
      "A",
      ~L,
      "id",
      "t",
      na.action = stats::na.exclude
    ),
    class = "survatr_bad_na_action"
  )
})

## Regression test for B2 (2026-04-22 critical review, round 1):
## `na.action = na.omit` (the default) combined with NA values in a
## predictor column silently misaligned `prep$X_fit` (post-NA-drop)
## against `fit_idx` (pre-NA-drop) in the sandwich IF chain, producing
## a subscript-out-of-bounds error. Fixed by rejecting NA in predictor
## columns upfront at surv_fit() via check_no_na_in_predictors(). Repro:
## `/tmp/survatr_repro_b2_na_drop.R`.
test_that("surv_fit rejects NA in predictor columns (B2)", {
  dt <- sim_constant_hazard(n = 200L, K = 4L, h = 0.1, seed = 441L)
  dt[, L := rnorm(.N)]
  dt_na_confounder <- data.table::copy(dt)
  dt_na_confounder[1L, L := NA]
  expect_error(
    surv_fit(dt_na_confounder, "Y", "A", ~L, "id", "t", time_formula = ~1),
    class = "survatr_na_in_predictors"
  )

  ## Also rejects NA in outcome, treatment, id, time.
  dt_na_outcome <- data.table::copy(dt)
  dt_na_outcome[1L, Y := NA]
  expect_error(
    surv_fit(dt_na_outcome, "Y", "A", ~L, "id", "t", time_formula = ~1),
    class = "survatr_na_in_predictors"
  )

  ## NA in censoring is NOT rejected (NA retains "uncensored" semantics).
  dt_na_cens <- data.table::copy(dt)
  dt_na_cens[, cens := 0L]
  dt_na_cens[1L, cens := NA]
  expect_silent(suppressWarnings(surv_fit(
    dt_na_cens,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "cens",
    time_formula = ~1
  )))
})

test_that("surv_fit rejects user-data collisions with reserved columns", {
  dt <- fixture_small_pp()
  dt[, .survatr_prev_event := 0L]
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t"),
    class = "survatr_reserved_col"
  )
})

test_that("print.survatr_fit emits a stable banner", {
  dt <- fixture_small_pp()
  fit <- suppressWarnings(surv_fit(
    dt,
    "Y",
    "A",
    ~L,
    "id",
    "t",
    censoring = "cens",
    time_formula = ~t
  ))
  expect_snapshot(print(fit))
})
