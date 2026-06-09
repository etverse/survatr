## Chunk 8: full matching rejection surface.
## survatr hard-rejects matching + survival on all entry routes:
##   1. estimator = "matching" or "match"       (Chunk 1 + Chunk 8)
##   2. method = "matching" / "match" in ...    (Chunk 8)
##   3. a `matchit` object passed as `data`     (Chunk 8, MatchIt-guarded)
##
## The correct tool for matched cohorts is survival::coxph with
## weights = match_weights and cluster = subclass directly on MatchIt
## output. survatr's pooled-logistic hazard on person-period data is
## structurally incompatible with matching weights (see hard-rules.md).

test_that("rejects estimator = 'matching' with classed error and redirect", {
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

test_that("rejects estimator = 'match' (alias) with classed error", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "match"),
    class = "survatr_matching_rejected"
  )
  expect_snapshot(
    surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "match"),
    error = TRUE
  )
})

test_that("rejects method = 'matching' in ... (causatr-style mis-call)", {
  ## A user might write surv_fit(..., method = "matching") by analogy with
  ## causatr's method= API. The check intercepts it before model dispatch.
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", method = "matching"),
    class = "survatr_matching_rejected"
  )
  expect_snapshot(
    surv_fit(dt, "Y", "A", ~L, "id", "t", method = "matching"),
    error = TRUE
  )
})

test_that("rejects method = 'match' in ... (causatr-style alias)", {
  dt <- fixture_small_pp()
  expect_error(
    surv_fit(dt, "Y", "A", ~L, "id", "t", method = "match"),
    class = "survatr_matching_rejected"
  )
  expect_snapshot(
    surv_fit(dt, "Y", "A", ~L, "id", "t", method = "match"),
    error = TRUE
  )
})

test_that("rejects a `matchit` object passed as data", {
  ## MatchIt's output object carries class "matchit"; survatr detects it
  ## before any column lookup so the error is informative rather than
  ## "column 'Y' not found".
  skip_if_not_installed("MatchIt")
  dt <- fixture_small_pp()
  ## Construct a real matchit object via a minimal 1:1 nearest-neighbor match.
  ## Muffle MatchIt's "fewer control units" balance warning; it is an artefact
  ## of the tiny fixture and irrelevant to the rejection being tested.
  mo <- withCallingHandlers(
    MatchIt::matchit(A ~ L, data = dt[!duplicated(dt$id)], method = "nearest"),
    warning = function(w) {
      if (grepl("control units", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_error(
    surv_fit(mo, "Y", "A", ~L, "id", "t"),
    class = "survatr_matching_rejected"
  )
  expect_snapshot(
    surv_fit(mo, "Y", "A", ~L, "id", "t"),
    error = TRUE
  )
})
