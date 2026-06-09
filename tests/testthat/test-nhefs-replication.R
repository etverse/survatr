## Track A acceptance test: Hernán & Robins Ch. 17 NHEFS replication.
## Reference: Hernán MA, Robins JM (2024). Causal Inference: What If. Ch. 17.
##   Targets (standardised survival curves, Program 17.2):
##     120-month survival ≈ 80.7% under qsmk = 1, ≈ 80.5% under qsmk = 0.
##     Risk difference ≈ 0.2% (essentially null), 95% CI spanning 0.
##
## The sandwich CI is compared against the published bootstrap CI as a ballpark
## (the two variance estimators need not match exactly).
##
## The nhefs_surv dataset was prepared by data-raw/nhefs_surv.R from
## causaldata::nhefs (rectangular 1629 × 120 person-period grid, 318 events).

## Lazy cache: the 195k-row GLM + sandwich computation is expensive; all tests
## in this file share a single fit computed on the first call.
.nhefs_cache <- new.env(parent = emptyenv())

.nhefs_load <- function(envir = .nhefs_cache) {
  if (exists("fit", envir = envir, inherits = FALSE)) {
    return(invisible(envir))
  }
  skip_on_cran()
  data("nhefs_surv", package = "survatr", envir = environment())
  nhefs_surv <- get("nhefs_surv", envir = environment())
  ## Pooled-logistic gcomp with the H&R Ch. 17 Program 17.2 hazard formula:
  ## quadratic time + treatment-by-time interactions.
  fit <- surv_fit(
    nhefs_surv,
    outcome = "event",
    treatment = "qsmk",
    confounders = ~ sex +
      race +
      age +
      I(age^2) +
      education +
      smokeintensity +
      I(smokeintensity^2) +
      smokeyrs +
      I(smokeyrs^2) +
      exercise +
      active +
      wt71 +
      I(wt71^2),
    id = "seqn",
    time = "time",
    time_formula = ~ time + I(time^2) + qsmk:time + qsmk:I(time^2),
    estimator = "gcomp"
  )
  res_surv <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(60L, 120L),
    type = "survival",
    ci_method = "sandwich"
  )
  res_rd <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = 120L,
    type = "risk_difference",
    ci_method = "sandwich"
  )
  assign("fit", fit, envir = envir)
  assign("res_surv", res_surv, envir = envir)
  assign("res_rd", res_rd, envir = envir)
  invisible(envir)
}

## Helper: extract s_hat at 120 months for one arm.
.s120 <- function(arm) {
  .nhefs_cache$res_surv$estimates[
    get("intervention") == arm & get("time") == 120L
  ]$s_hat
}

test_that("120-month survival under qsmk=1 ≈ 80.7% (H&R Ch.17 target)", {
  skip_on_cran()
  .nhefs_load()
  expect_equal(.s120("a1"), 0.807, tolerance = 0.03)
})

test_that("120-month survival under qsmk=0 ≈ 80.5% (H&R Ch.17 target)", {
  skip_on_cran()
  .nhefs_load()
  expect_equal(.s120("a0"), 0.805, tolerance = 0.03)
})

test_that("risk difference at 120 months ≈ 0.2% (essentially null)", {
  skip_on_cran()
  .nhefs_load()
  rd_row <- .nhefs_cache$res_rd$contrasts[get("time") == 120L]
  expect_equal(rd_row$estimate, 0.002, tolerance = 0.01)
  ## Sandwich CI must span zero (null treatment effect).
  expect_lt(rd_row$ci_lower, 0)
  expect_gt(rd_row$ci_upper, 0)
})

test_that("sandwich CI width is in the ballpark of H&R bootstrap CI (≈ 7.8 pp)", {
  skip_on_cran()
  .nhefs_load()
  rd_row <- .nhefs_cache$res_rd$contrasts[get("time") == 120L]
  ci_width <- rd_row$ci_upper - rd_row$ci_lower
  ## H&R published bootstrap CI: −4.1% to 3.7% = 7.8 pp wide. Guard against
  ## implausible widths (< 4 pp or > 20 pp) that would indicate a variance bug.
  expect_gt(ci_width, 0.04)
  expect_lt(ci_width, 0.20)
})

test_that("unadjusted KM at 120 months is in the same ballpark as standardised curves", {
  skip_on_cran()
  skip_if_not_installed("survival")
  data("nhefs_surv", package = "survatr", envir = environment())
  nhefs_surv <- get("nhefs_surv", envir = environment())
  ## Collapse rectangular PP to one-row-per-id for KM.
  per_id <- nhefs_surv[,
    {
      ev_row <- which(event == 1L)
      if (length(ev_row)) {
        list(t = time[ev_row[1L]], status = 1L)
      } else {
        list(t = 120L, status = 0L)
      }
    },
    by = "seqn"
  ]
  km <- survival::survfit(survival::Surv(t, status) ~ 1, data = per_id)
  km_120 <- summary(km, times = 120L, extend = TRUE)$surv
  ## Unadjusted survival should be in the same 75–90% window as the standardised
  ## curves (roughly null treatment effect, so little adjustment).
  expect_gt(km_120, 0.75)
  expect_lt(km_120, 0.90)
})
