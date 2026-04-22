## Validation / rejection-path regression tests.

make_fit <- function() {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 101L)
  surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
}

test_that("contrast rejects empty or unnamed interventions", {
  fit <- make_fit()
  expect_error(
    contrast(fit, interventions = list(), times = 1:5, type = "survival"),
    class = "survatr_bad_interventions"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(causatr::static(0)),
      times = 1:5,
      type = "survival"
    ),
    class = "survatr_bad_interventions"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a = causatr::static(0), a = causatr::static(1)),
      times = 1:5,
      type = "survival"
    ),
    class = "survatr_bad_interventions"
  )
})

test_that("contrast rejects non-intervention list elements", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = 0, a1 = 1),
      times = 1:5,
      type = "survival"
    ),
    class = "survatr_bad_interventions"
  )
})

## Regression test for R3 (2026-04-22 critical review, round 1):
## `validate_times` previously used `is.numeric(times)`, which returns
## FALSE for Date / POSIXct / difftime. Users whose hazard model uses a
## real-world timestamp column could not evaluate contrasts. Relaxed
## to accept numeric and the three common time-like classes, delegating
## set-membership to `setdiff()` which works across them.
test_that("contrast accepts Date / POSIXct times when fit uses them (R3)", {
  dt <- sim_constant_hazard(n = 100L, K = 3L, h = 0.1, seed = 481L)
  dt[, t := as.Date("2020-01-01") + (t - 1L) * 30L]
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
  expect_silent(contrast(
    fit,
    list(a0 = causatr::static(0)),
    times = fit$time_grid,
    type = "survival"
  ))
})

test_that("contrast rejects bad `times`", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = numeric(0),
      type = "survival"
    ),
    class = "survatr_bad_times"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = c(1, NA),
      type = "survival"
    ),
    class = "survatr_bad_times"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = "1",
      type = "survival"
    ),
    class = "survatr_bad_times"
  )
})

test_that("contrast rejects extrapolation beyond observed time grid", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = c(1, 2, 99),
      type = "survival"
    ),
    class = "survatr_time_extrapolation"
  )
})

test_that("contrast rejects unknown `type`", {
  fit <- make_fit()
  ## `match.arg` raises a simpleError, not a classed one -- regression-test
  ## that the message mentions the offending value.
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      type = "bogus"
    ),
    "should be one of"
  )
})

test_that("contrast rejects a `reference` that is not in interventions", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 1:5,
      type = "risk_difference",
      reference = "missing"
    ),
    class = "survatr_bad_reference"
  )
})

test_that("contrast rejects bad conf_level", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      type = "survival",
      ci_method = "sandwich",
      conf_level = 1.5
    ),
    class = "survatr_bad_conf_level"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      type = "survival",
      ci_method = "sandwich",
      conf_level = 0
    ),
    class = "survatr_bad_conf_level"
  )
})

test_that("contrast rejects unknown ci_method", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      type = "survival",
      ci_method = "bogus"
    ),
    class = "survatr_bad_ci_method"
  )
})
