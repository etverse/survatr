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
  wide <- data.table::data.table(
    id = 1:5,
    t = 1L,
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

test_that("prepare_pp_data rejects duplicated (id, time) rows", {
  dup <- rbind(fixture_small_pp(), fixture_small_pp()[1L])
  expect_error(
    prepare_pp_data(dup, "Y", "A", "id", "t"),
    class = "survatr_duplicate_pp_row"
  )
})
