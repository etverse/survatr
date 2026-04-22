make_res <- function(type = "risk_difference", ci = "sandwich") {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 401L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 1:5,
    type = type,
    reference = "a0",
    ci_method = ci
  )
}

test_that("tidy returns a long data.frame with expected columns", {
  res <- make_res()
  out <- tidy(res)
  expect_s3_class(out, "data.frame")
  expect_setequal(
    names(out),
    c(
      "intervention",
      "contrast",
      "time",
      "estimand",
      "estimate",
      "se",
      "ci_lower",
      "ci_upper"
    )
  )
  ## 2 interventions x 5 times + 1 contrast x 5 times = 15 rows.
  expect_equal(nrow(out), 15L)
})

test_that("tidy which = 'estimates' drops contrast rows", {
  res <- make_res()
  out <- tidy(res, which = "estimates")
  expect_true(all(is.na(out$contrast)))
  expect_equal(nrow(out), 10L) ## 2 interventions x 5 times
})

test_that("tidy which = 'contrasts' drops intervention rows", {
  res <- make_res()
  out <- tidy(res, which = "contrasts")
  expect_true(all(is.na(out$intervention)))
  expect_equal(nrow(out), 5L)
})

test_that("tidy conf.int = FALSE drops CI columns", {
  res <- make_res()
  out <- tidy(res, conf.int = FALSE)
  expect_false("ci_lower" %in% names(out))
  expect_false("ci_upper" %in% names(out))
})

test_that("tidy on a curve-only type returns empty contrast rows", {
  res <- make_res(type = "survival", ci = "sandwich")
  out <- tidy(res)
  ## `type = survival` has empty contrasts table; `which = all` reduces to
  ## estimates only.
  expect_equal(nrow(out), 10L)
  expect_true(all(is.na(out$contrast)))
})

test_that("tidy rejects invalid conf.int argument", {
  res <- make_res()
  expect_error(tidy(res, conf.int = NA), class = "survatr_bad_conf_int")
  expect_error(tidy(res, conf.int = "yes"), class = "survatr_bad_conf_int")
})
