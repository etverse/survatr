## Competing risks: cause-specific hazards + CIF point estimates, the
## F1 + F2 + S = 1 identity, cause subsetting, and the rejection paths.
## Oracles in helper-cr-oracle.R: closed-form constant-hazard CIF + the
## Aalen-Johansen estimator (survival::survfit).

test_that("competing-risks fit has the expected structure", {
  pp <- sim_two_cause_constant_hazard(n = 400L, K = 5L, seed = 3L)
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~1
  )
  expect_s3_class(fit, "survatr_fit")
  expect_null(fit$model)
  expect_identical(fit$causes, c(1L, 2L))
  expect_named(fit$cause_models, c("1", "2"))
  expect_s3_class(fit$cause_models[["1"]], "glm")
  expect_identical(fit$track, "A")
  expect_identical(fit$competing, "event")
  ## Both cause models fit on the SAME shared all-cause risk set.
  expect_identical(
    stats::nobs(fit$cause_models[["1"]]),
    stats::nobs(fit$cause_models[["2"]])
  )
})

test_that("CIF point estimates match the closed-form two-cause oracle", {
  ## Constant hazards, no treatment / covariate effect: F_j(t) =
  ## (h_j / (h1 + h2)) * (1 - (1 - h1 - h2)^t) exactly (up to MC error).
  h1 <- 0.06
  h2 <- 0.04
  pp <- sim_two_cause_constant_hazard(
    n = 6000L,
    K = 6L,
    h1 = h1,
    h2 = h2,
    seed = 5L
  )
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~1
  )
  times <- c(2L, 4L, 6L)
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = times,
    type = "cif"
  )
  truth <- analytic_cr(times, h1, h2)
  f1 <- res$estimates[get("cause") == 1L][order(time)]
  f2 <- res$estimates[get("cause") == 2L][order(time)]
  expect_equal(f1$cif_hat, truth$F1, tolerance = 0.02)
  expect_equal(f2$cif_hat, truth$F2, tolerance = 0.02)

  ## All-cause survival from the summed hazards matches (1 - h1 - h2)^t.
  res_s <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = times,
    type = "survival"
  )
  expect_equal(res_s$estimates[order(time)]$s_hat, truth$S, tolerance = 0.02)
})

test_that("the F1 + F2 + S = 1 identity holds numerically", {
  pp <- sim_two_cause_constant_hazard(n = 1500L, K = 5L, seed = 9L)
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~ factor(t)
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  times <- 1:5
  cif <- contrast(fit, interventions = ivs, times = times, type = "cif")
  surv <- contrast(fit, interventions = ivs, times = times, type = "survival")
  ## Sum the per-cause CIFs within (intervention, time), add S, expect 1.
  fsum <- cif$estimates[,
    list(fsum = sum(cif_hat)),
    by = c("intervention", "time")
  ]
  m <- merge(
    fsum,
    surv$estimates[, c("intervention", "time", "s_hat")],
    by = c("intervention", "time")
  )
  expect_equal(m$fsum + m$s_hat, rep(1, nrow(m)), tolerance = 1e-12)
})

test_that("gcomp CIF matches the Aalen-Johansen estimator (saturated)", {
  skip_if_not_installed("survival")
  ## A saturated baseline hazard (~ factor(t)) reproduces the per-period
  ## empirical cause-specific hazards AJ uses, so on a no-confounding DGP the
  ## standardized gcomp CIF and the AJ CIF agree closely.
  pp <- sim_two_cause_constant_hazard(n = 4000L, K = 5L, seed = 12L)
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~ factor(t)
  )
  times <- 1:5
  res <- contrast(
    fit,
    interventions = list(a0 = causatr::static(0)),
    times = times,
    type = "cif"
  )
  aj <- data.table::as.data.table(aj_cif_oracle(pp, times))
  ## Merge gcomp CIF and AJ CIF on (cause, time) -- robust to row ordering.
  cmp <- merge(
    res$estimates[, c("cause", "time", "cif_hat")],
    aj,
    by = c("cause", "time")
  )
  expect_equal(nrow(cmp), length(times) * 2L)
  expect_equal(cmp$cif_hat, cmp$cif, tolerance = 0.02)
})

test_that("cause = subsets the reported causes", {
  pp <- sim_two_cause_constant_hazard(n = 600L, K = 4L, seed = 4L)
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~1
  )
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(2L, 4L),
    type = "cif_difference",
    cause = 1L
  )
  expect_identical(sort(unique(res$estimates$cause)), 1L)
  expect_identical(sort(unique(res$contrasts$cause)), 1L)
})

test_that("competing-risks bootstrap populates CIs around the point", {
  pp <- sim_two_cause_constant_hazard(
    n = 500L,
    K = 4L,
    seed = 6L,
    beta1_A = -0.5,
    gamma = 0.4
  )
  fit <- surv_fit(
    pp,
    outcome = "event",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~ factor(t)
  )
  res <- suppressMessages(contrast(
    fit,
    interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
    times = c(2L, 4L),
    type = "cif_difference",
    cause = 1L,
    ci_method = "bootstrap",
    n_boot = 80L,
    seed = 1L
  ))
  expect_false(anyNA(res$contrasts$se))
  expect_true(all(res$contrasts$ci_lower <= res$contrasts$estimate))
  expect_true(all(res$contrasts$estimate <= res$contrasts$ci_upper))
})

test_that("competing-risks rejection paths fire with classed errors", {
  pp <- sim_two_cause_constant_hazard(n = 200L, K = 3L, seed = 2L)
  base <- list(
    data = pp,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    time_formula = ~1
  )
  ## IPW / ICE + competing => estimator rejection.
  expect_error(
    do.call(surv_fit, c(base, list(competing = "event", estimator = "ipw"))),
    class = "survatr_competing_estimator"
  )
  ## outcome != competing => misuse.
  expect_error(
    surv_fit(
      pp,
      outcome = "A",
      treatment = "A",
      confounders = ~1,
      id = "id",
      time = "t",
      competing = "event"
    ),
    class = "survatr_competing_misuse"
  )
  ## Single positive cause => misuse (not a competing-risks problem).
  pp1 <- data.table::copy(pp)
  pp1[event == 2L, event := 1L]
  expect_error(
    do.call(
      surv_fit,
      c(
        list(data = pp1),
        base[-1L],
        list(competing = "event")
      )
    ),
    class = "survatr_competing_misuse"
  )
  ## Negative / non-integer cause code => bad competing column.
  pp_bad <- data.table::copy(pp)
  pp_bad[1L, event := -1L]
  expect_error(
    do.call(
      surv_fit,
      c(
        list(data = pp_bad),
        base[-1L],
        list(competing = "event")
      )
    ),
    class = "survatr_bad_competing"
  )
})

test_that("CIF estimands require a competing-risks fit and vice versa", {
  pp <- sim_two_cause_constant_hazard(n = 300L, K = 3L, seed = 7L)
  ## Single-event fit asked for a CIF estimand.
  se_fit <- surv_fit(
    pp[event == 2L, event := 0L][],
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    time_formula = ~1
  )
  expect_error(
    contrast(
      se_fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 2L,
      type = "cif_difference"
    ),
    class = "survatr_competing_type"
  )
  ## Competing fit asked for a single-event-only contrast (rmst).
  pp2 <- sim_two_cause_constant_hazard(n = 300L, K = 3L, seed = 7L)
  cr_fit <- surv_fit(
    pp2,
    outcome = "event",
    treatment = "A",
    confounders = ~1,
    id = "id",
    time = "t",
    competing = "event",
    time_formula = ~1
  )
  expect_error(
    contrast(
      cr_fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 2L,
      type = "rmst_difference"
    ),
    class = "survatr_competing_type"
  )
  ## Unknown cause label.
  expect_error(
    contrast(
      cr_fit,
      interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
      times = 2L,
      type = "cif_difference",
      cause = 9L
    ),
    class = "survatr_bad_cause"
  )
})
