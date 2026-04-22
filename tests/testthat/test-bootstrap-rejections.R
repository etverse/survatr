make_fit_small <- function() {
  dt <- sim_constant_hazard(n = 200L, K = 4L, h = 0.1, seed = 301L)
  surv_fit(dt, "Y", "A", ~1, "id", "t", time_formula = ~1)
}

test_that("contrast rejects single-intervention with a pairwise-contrast type", {
  fit <- make_fit_small()
  for (ty in c("risk_difference", "risk_ratio", "rmst_difference")) {
    expect_error(
      contrast(
        fit,
        interventions = list(a0 = causatr::static(0)),
        times = 1:4,
        type = ty
      ),
      class = "survatr_bad_interventions"
    )
  }
})

test_that("contrast rejects bad n_boot", {
  fit <- make_fit_small()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      n_boot = 0L
    ),
    class = "survatr_bad_n_boot"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      n_boot = 10.5
    ),
    class = "survatr_bad_n_boot"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      n_boot = -10L
    ),
    class = "survatr_bad_n_boot"
  )
})

test_that("contrast rejects bad boot_ci", {
  fit <- make_fit_small()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      boot_ci = "bca"
    ),
    class = "survatr_bad_boot_ci"
  )
})

test_that("contrast rejects bad parallel / ncpus", {
  fit <- make_fit_small()
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      parallel = "bogus"
    ),
    class = "survatr_bad_parallel"
  )
  expect_error(
    contrast(
      fit,
      interventions = list(a0 = causatr::static(0)),
      times = 1:4,
      type = "survival",
      ci_method = "bootstrap",
      parallel = "no",
      ncpus = 0L
    ),
    class = "survatr_bad_parallel"
  )
})
