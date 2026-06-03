test_that("forrest() runs on a contrast-shaped result", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 601L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference",
    reference = "a0",
    ci_method = "sandwich"
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(
    {
      grDevices::dev.off()
      unlink(tmp)
    },
    add = TRUE
  )
  expect_invisible(forrest(res, t_ref = 3))
})

test_that("forrest() rejects a t_ref outside time_grid", {
  dt <- sim_constant_hazard(n = 300L, K = 5L, h = 0.1, seed = 603L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference",
    reference = "a0",
    ci_method = "sandwich"
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(
    {
      grDevices::dev.off()
      unlink(tmp)
    },
    add = TRUE
  )
  expect_error(forrest(res, t_ref = 99), class = "survatr_bad_t_ref")
  expect_error(forrest(res, t_ref = c(1, 2)), class = "survatr_bad_t_ref")
})

test_that("forrest() rejects a curve-only result", {
  dt <- sim_constant_hazard(n = 300L, K = 5L, h = 0.1, seed = 607L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:5,
    type = "survival",
    ci_method = "sandwich"
  )
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(
    {
      grDevices::dev.off()
      unlink(tmp)
    },
    add = TRUE
  )
  expect_error(
    forrest(res, t_ref = 3),
    class = "survatr_forrest_wrong_type"
  )
})

test_that("forrest() distinguishes off-grid t_ref from no-contrasts-at-t_ref", {
  dt <- sim_constant_hazard(n = 300L, K = 5L, h = 0.1, seed = 609L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = "risk_difference",
    reference = "a0",
    ci_method = "sandwich"
  )
  # t_ref = 3 stays in time_grid (so the off-grid guard passes), but the
  # contrasts table has no rows there -- the empty-contrasts case a
  # single-intervention result would produce. It must abort with a class
  # distinct from survatr_bad_t_ref so callers can tell the two apart.
  res$contrasts <- res$contrasts[get("time") != 3]
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(
    {
      grDevices::dev.off()
      unlink(tmp)
    },
    add = TRUE
  )
  expect_error(
    forrest(res, t_ref = 3),
    class = "survatr_forrest_no_contrasts"
  )
})
