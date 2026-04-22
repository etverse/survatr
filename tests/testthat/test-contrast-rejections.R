## Validation / rejection-path regression tests.

make_fit <- function() {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 101L)
  surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
}

test_that("contrast rejects empty or unnamed interventions", {
  fit <- make_fit()
  expect_error(
    contrast(fit, interventions = list(), times = 1:5),
    class = "survatr_bad_interventions"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(causatr::static(0)),
      times = 1:5
    ),
    class = "survatr_bad_interventions"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a = causatr::static(0), a = causatr::static(1)),
      times = 1:5
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
      times = 1:5
    ),
    class = "survatr_bad_interventions"
  )
})

test_that("contrast rejects bad `times`", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = numeric(0)
    ),
    class = "survatr_bad_times"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = c(1, NA)
    ),
    class = "survatr_bad_times"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = "1"
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
      times = c(1, 2, 99)
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

test_that("contrast rejects ci_method in {sandwich, bootstrap} with a pointer", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 1:5,
      ci_method = "sandwich"
    ),
    class = "survatr_ci_not_available"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 1:5,
      ci_method = "bootstrap"
    ),
    class = "survatr_ci_not_available"
  )
})

test_that("contrast rejects unknown ci_method", {
  fit <- make_fit()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      ci_method = "bogus"
    ),
    class = "survatr_bad_ci_method"
  )
})
