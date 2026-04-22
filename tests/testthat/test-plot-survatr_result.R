## plot() tests: we are not using vdiffr so these are smoke-only --
## assert that the call runs without error under each `type` and returns
## the object invisibly. Graphical regressions are caught by the
## type-level column coverage in other tests.

skip_plot <- function() {
  ## Guard against non-interactive graphics failures on CI. Open a PDF
  ## device so no X11 / Quartz is needed.
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(file = tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  on.exit(unlink(tmp), add = TRUE)
  invisible(TRUE)
}

test_that("plot() runs on survival / risk / rmst curves", {
  skip_plot()
  dt <- sim_constant_hazard(n = 300L, K = 5L, h = 0.1, seed = 501L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  for (type in c("survival", "risk", "rmst")) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(file = tmp)
    res <- contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:5,
      type = type,
      ci_method = "sandwich"
    )
    expect_invisible(plot(res))
    grDevices::dev.off()
    unlink(tmp)
  }
})

test_that("plot() runs on risk_difference / risk_ratio / rmst_difference", {
  dt <- sim_constant_hazard(n = 500L, K = 5L, h = 0.1, seed = 503L)
  fit <- surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
  for (type in c("risk_difference", "risk_ratio", "rmst_difference")) {
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(file = tmp)
    res <- contrast(
      fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 1:5,
      type = type,
      reference = "a0",
      ci_method = "sandwich"
    )
    expect_invisible(plot(res))
    grDevices::dev.off()
    unlink(tmp)
  }
})

test_that("plot() which = 'contrasts' on a curve-only result aborts cleanly", {
  dt <- sim_constant_hazard(n = 300L, K = 5L, h = 0.1, seed = 507L)
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
    plot(res, which = "contrasts"),
    class = "survatr_plot_no_contrasts"
  )
})
