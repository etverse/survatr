## Cause-aware S3 methods (print / tidy / plot / forrest) on competing-risks
## results. The single-event S3 tests live in test-{print,tidy,plot,forrest}-*;
## these pin the cause dimension the competing-risks path adds.

# Small shared CR fit (cheap; built once per test for isolation).
cr_fit <- function(n = 500L, seed = 5L) {
  pp <- sim_two_cause_constant_hazard(
    n = n,
    K = 4L,
    seed = seed,
    beta1_A = -0.5,
    gamma = 0.4
  )
  surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~ factor(t)
  )
}

cr_result <- function(type, cause = NULL) {
  suppressMessages(contrast(
    cr_fit(),
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(2L, 4L),
    type = type,
    cause = cause,
    reference = "a1",
    ci_method = "sandwich"
  ))
}

test_that("print shows the truncation-by-death caveat for CIF contrasts", {
  res <- cr_result("cif_difference", cause = 1L)
  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(out, "truncation by death")
  expect_match(out, "cif_difference")
})

test_that("print omits the caveat for all-cause survival on a CR fit", {
  res <- cr_result("survival")
  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_false(grepl("truncation by death", out))
})

test_that("tidy carries the cause column for CIF results", {
  res <- cr_result("cif_difference", cause = c(1L, 2L))
  td <- tidy(res)
  expect_true("cause" %in% names(td))
  expect_setequal(stats::na.omit(unique(td$cause)), c(1L, 2L))
  ## Intervention rows use cif_hat; contrast rows use the contrast type.
  expect_true("cif_hat" %in% td$estimand)
  expect_true("cif_difference" %in% td$estimand)
})

test_that("tidy of all-cause survival on a CR fit keeps cause = NA", {
  res <- cr_result("survival")
  td <- tidy(res)
  ## The estimates carry a cause column (all NA for all-cause rows).
  expect_true("cause" %in% names(td))
  expect_true(all(is.na(td$cause)))
})

test_that("plot renders CIF curves and CIF-difference contrasts", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(cr_result("cif", cause = c(1L, 2L))))
  expect_invisible(plot(cr_result("cif_difference", cause = c(1L, 2L))))
  expect_invisible(plot(cr_result("survival")))
})

test_that("forrest renders one row per (contrast, cause) at t_ref", {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- cr_result("cif_ratio", cause = c(1L, 2L))
  expect_invisible(forrest(res, t_ref = 4L))
  ## A curve-only CIF result has no contrasts -> forrest must reject.
  expect_error(
    forrest(cr_result("cif", cause = 1L), t_ref = 4L),
    class = "survatr_forrest_wrong_type"
  )
})
