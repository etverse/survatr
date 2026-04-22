test_that("prepare_pp_data coerces data.frame to data.table", {
  df <- as.data.frame(fixture_small_pp())
  out <- prepare_pp_data(df, "Y", "A", "id", "t")
  expect_s3_class(out, "data.table")
  expect_equal(nrow(out), 15L)
})

test_that("prepare_pp_data copies input data.table (no in-place mutation)", {
  dt <- fixture_small_pp()
  snap <- data.table::copy(dt)
  out <- prepare_pp_data(dt, "Y", "A", "id", "t")
  ## Mutate the returned object; original should be unchanged.
  out[, .survatr_marker := 1L]
  expect_identical(names(dt), names(snap))
  expect_false(".survatr_marker" %in% names(dt))
})

test_that("prepare_pp_data sorts by (id, time)", {
  dt <- fixture_small_pp()
  shuffled <- dt[sample(.N)]
  out <- prepare_pp_data(shuffled, "Y", "A", "id", "t")
  expect_equal(out$id, rep(1:5, each = 3L))
  expect_equal(out$t, rep(1:3, times = 5L))
})

test_that("prepare_pp_data rejects missing columns", {
  dt <- fixture_small_pp()
  expect_error(
    prepare_pp_data(dt, "Ymissing", "A", "id", "t"),
    class = "survatr_col_not_found"
  )
  expect_snapshot(
    prepare_pp_data(dt, "Y", "A", "id", "t", censoring = "cens_missing"),
    error = TRUE
  )
})

test_that("prepare_pp_data rejects wide (one row per id) input", {
  ## A wide fixture: 5 ids, each with one row, but the time grid has
  ## n_times > 1 so at least some id is missing most time values. The
  ## message distinguishes "wide" from "ragged" at the n_times == 1
  ## degenerate edge case.
  wide <- data.table::data.table(
    id = c(1, 2, 3, 4, 5),
    t = c(1, 2, 3, 1, 2),
    A = c(1, 0, 1, 0, 1),
    Y = c(0, 1, 0, 1, 0)
  )
  expect_error(
    prepare_pp_data(wide, "Y", "A", "id", "t"),
    class = "survatr_not_person_period"
  )
  expect_snapshot(
    prepare_pp_data(wide, "Y", "A", "id", "t"),
    error = TRUE
  )
})

## Regression test for B3 / S3 (2026-04-22 critical review, round 1):
## Ragged PP -- ids missing a row at some time in the unique-time grid --
## previously crashed the sandwich IF chain with a subscript-OOB at
## compute_survival_if_matrix line 111. Now rejected at prepare_pp_data
## with a dedicated class `survatr_ragged_pp`. The classic "one id had
## the event at t=1 and post-event rows were dropped" case is the target
## scenario. Repro: `/tmp/survatr_repro_b3_ragged_pp.R`.
test_that("prepare_pp_data rejects ragged PP (B3 / S3)", {
  ## Standard ragged case: id 1 event at t=1, only (1, 1) present.
  rect <- data.table::CJ(id = 1:5, t = 1:3)
  rect[, A := 1L]
  rect[, Y := 0L]
  data.table::setcolorder(rect, c("id", "t", "A", "Y"))
  ragged <- rect[-2L] ## drop (1, 2), leaving 14 rows instead of 15
  expect_error(
    prepare_pp_data(ragged, "Y", "A", "id", "t"),
    class = "survatr_ragged_pp"
  )

  ## Single-row ids in a 1-period study (n_times = 1) are trivially
  ## rectangular and must be accepted. This replaces the old hard
  ## "one-row-per-id = wide" heuristic which conflated ragged and wide.
  one_period <- data.table::data.table(id = 1:5, t = 1L, A = 0L, Y = 0L)
  expect_silent(prepare_pp_data(one_period, "Y", "A", "id", "t"))
})

test_that("prepare_pp_data rejects duplicated (id, time) rows", {
  dup <- rbind(fixture_small_pp(), fixture_small_pp()[1L])
  expect_error(
    prepare_pp_data(dup, "Y", "A", "id", "t"),
    class = "survatr_duplicate_pp_row"
  )
})
