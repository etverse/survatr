test_that("is_uncensored treats NA and 0 as at-risk, anything else as censored", {
  dt <- data.table::data.table(
    cens = c(0L, 1L, NA, 0L, 2L)
  )
  expect_equal(is_uncensored(dt, "cens"), c(TRUE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(is_uncensored(dt, NULL), rep(TRUE, 5L))
})

test_that("build_risk_set drops rows at/after the first event per id", {
  dt <- prepare_pp_data(fixture_small_pp(), "Y", "A", "id", "t")
  fit_rows <- build_risk_set(dt, outcome = "Y", id = "id", censoring = NULL)

  ## Expected mask (censoring ignored in this call):
  ## id 1: TRUE, TRUE, FALSE (event at t=2; t=3 at/after first event)
  ## id 2: TRUE, TRUE, TRUE (no event)
  ## id 3: TRUE, TRUE, TRUE (cens ignored; no event in Y)
  ## id 4: TRUE, FALSE, FALSE (event at t=1)
  ## id 5: TRUE, TRUE, TRUE
  expect_equal(
    fit_rows,
    c(
      TRUE,
      TRUE,
      FALSE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      FALSE,
      FALSE,
      TRUE,
      TRUE,
      TRUE
    )
  )
  expect_equal(sum(fit_rows), 12L)
})

test_that("build_risk_set drops rows at/after the first censor when censoring is passed", {
  dt <- prepare_pp_data(
    fixture_small_pp(),
    "Y",
    "A",
    "id",
    "t",
    censoring = "cens"
  )
  fit_rows <- build_risk_set(dt, outcome = "Y", id = "id", censoring = "cens")

  ## id 3 censored at t=2: at-risk row (3, t=1) only.
  ## All others as in the no-censor case except id 3.
  expect_equal(
    fit_rows,
    c(
      TRUE,
      TRUE,
      FALSE,
      TRUE,
      TRUE,
      TRUE,
      TRUE,
      FALSE,
      FALSE, ## id 3: t=2 censored; t=3 dropped by prev_cens
      TRUE,
      FALSE,
      FALSE,
      TRUE,
      TRUE,
      TRUE
    )
  )
  expect_equal(sum(fit_rows), 10L)
})

test_that("build_risk_set writes the internal columns onto data (mutation by design)", {
  dt <- prepare_pp_data(
    fixture_small_pp(),
    "Y",
    "A",
    "id",
    "t",
    censoring = "cens"
  )
  build_risk_set(dt, outcome = "Y", id = "id", censoring = "cens")
  expect_true(".survatr_prev_event" %in% names(dt))
  expect_true(".survatr_prev_cens" %in% names(dt))
})
