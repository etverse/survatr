test_that("check_weights accepts NULL and zero weights", {
  expect_silent(check_weights(NULL, 10L))
  expect_silent(check_weights(rep(0, 10L), 10L))
  expect_silent(check_weights(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), 10L))
})

test_that("check_weights rejects non-numeric / mis-sized / NA / Inf / negative", {
  expect_snapshot(check_weights(c("a", "b", "c"), 3L), error = TRUE)
  expect_snapshot(check_weights(c(1, 2, 3), 5L), error = TRUE)
  expect_snapshot(check_weights(c(1, NA, 3), 3L), error = TRUE)
  expect_snapshot(check_weights(c(1, Inf, 3), 3L), error = TRUE)
  expect_snapshot(check_weights(c(1, NaN, 3), 3L), error = TRUE)
  expect_snapshot(check_weights(c(1, -2, 3), 3L), error = TRUE)
})

test_that("check_weights errors carry the survatr_bad_weights class", {
  expect_error(check_weights(c(1, -2, 3), 3L), class = "survatr_bad_weights")
  expect_error(check_weights(c(1, 2), 3L), class = "survatr_bad_weights")
})

test_that("check_dots_na_action accepts na.omit / na.fail / absent", {
  expect_silent(check_dots_na_action())
  expect_silent(check_dots_na_action(na.action = stats::na.omit))
  expect_silent(check_dots_na_action(na.action = stats::na.fail))
  expect_silent(check_dots_na_action(na.action = "na.omit"))
  expect_silent(check_dots_na_action(na.action = "na.fail"))
})

test_that("check_dots_na_action rejects na.exclude (function and string)", {
  expect_error(
    check_dots_na_action(na.action = stats::na.exclude),
    class = "survatr_bad_na_action"
  )
  expect_error(
    check_dots_na_action(na.action = "na.exclude"),
    class = "survatr_bad_na_action"
  )
  expect_snapshot(
    check_dots_na_action(na.action = stats::na.exclude),
    error = TRUE
  )
})

## Regression test for S1 (2026-04-22 critical review, round 1):
## `.cf_hazard` / `.cf_surv` are now also reserved -- used as temporary
## columns on counterfactual PP copies in `compute_survival_curve()` and
## `compute_survival_if_matrix()`. A user data column with either name
## would have been silently overwritten in the internal copy; the
## reserved-col guard at the boundary now catches that upfront.
test_that("check_reserved_cols rejects .cf_hazard / .cf_surv (S1)", {
  ok <- data.table::data.table(id = 1L, t = 1L, y = 0L)
  bad_cf_hazard <- data.table::copy(ok)
  bad_cf_hazard[, .cf_hazard := 0]
  expect_error(
    check_reserved_cols(bad_cf_hazard),
    class = "survatr_reserved_col"
  )

  bad_cf_surv <- data.table::copy(ok)
  bad_cf_surv[, .cf_surv := 0]
  expect_error(
    check_reserved_cols(bad_cf_surv),
    class = "survatr_reserved_col"
  )
})

test_that("check_reserved_cols rejects .survatr_prev_event / .survatr_prev_cens", {
  ok <- data.table::data.table(id = 1L, t = 1L, y = 0L)
  expect_silent(check_reserved_cols(ok))

  bad1 <- data.table::copy(ok)
  bad1[, .survatr_prev_event := 0L]
  expect_error(check_reserved_cols(bad1), class = "survatr_reserved_col")

  bad2 <- data.table::copy(ok)
  bad2[, .survatr_prev_cens := 0L]
  expect_error(check_reserved_cols(bad2), class = "survatr_reserved_col")

  bad_both <- data.table::copy(ok)
  bad_both[, `:=`(.survatr_prev_event = 0L, .survatr_prev_cens = 0L)]
  expect_snapshot(check_reserved_cols(bad_both), error = TRUE)
})
