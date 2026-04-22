test_that("type = 'survival' returns s_hat, empty contrasts, no reference", {
  dt <- sim_constant_hazard(n = 1000L, K = 6L, h = 0.08, seed = 41L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:6,
    type = "survival"
  )
  expect_true("s_hat" %in% names(res$estimates))
  expect_equal(nrow(res$estimates), 6L)
  expect_equal(nrow(res$contrasts), 0L)
  expect_null(res$reference)
})

test_that("type = 'risk' returns risk_hat and empty contrasts", {
  dt <- sim_constant_hazard(n = 1000L, K = 6L, h = 0.08, seed = 43L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:6,
    type = "risk"
  )
  expect_true("risk_hat" %in% names(res$estimates))
  expect_equal(nrow(res$contrasts), 0L)
  expect_equal(res$type, "risk")
  expect_null(res$reference)
})

test_that("survival curve under static(0) matches (1 - h)^t end-to-end", {
  h <- 0.05
  dt <- sim_constant_hazard(n = 5000L, K = 10L, h = h, seed = 47L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = 1:10,
    type = "survival"
  )
  expect_equal(res$estimates$s_hat, (1 - h)^(1:10), tolerance = 0.01)
})
