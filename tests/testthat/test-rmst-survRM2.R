## KM-based RMST cross-check: survatr gcomp RMST vs survRM2::rmst2().
##
## Oracle rationale:
##   survRM2 computes the restricted mean survival time via the area under the
##   KM curve (continuous step function) while survatr computes the
##   trapezoidal integral of the pooled-logistic survival curve. For a
##   constant-hazard DGP with no confounders and small per-period hazard
##   (h = 0.05), the two agree closely: the pooled-logistic model recovers
##   the discrete (1-h)^t curve, whose trapezoidal integral approximates the
##   KM area with a bias of order h²/12 per interval — well below the 0.05
##   tolerance used here. This is a sanity check that the RMST *scale* is
##   correct; tight precision pins live in test-rmst.R (closed-form
##   trapezoidal oracle) and test-gcomp-delicatessen.R (delicatessen).
##
## DGP:
##   n = 3000, K = 20, h = 0.05 (constant, no treatment effect, no covariates).
##   Both arms expected to have the same true RMST (no treatment effect; gcomp
##   with confounders = ~1 and a 50/50 treatment gives ~identical arm curves).

test_that("gcomp RMST agrees with survRM2 (unadjusted, constant hazard)", {
  skip_if_not_installed("survRM2")
  skip_on_cran()

  dat <- sim_constant_hazard(n = 3000L, K = 20L, h = 0.05, seed = 42L)

  ## Constant hazard: time_formula = ~1 matches the DGP exactly.
  fit <- surv_fit(
    dat,
    "Y",
    "A",
    ~1,
    "id",
    "t",
    time_formula = ~1,
    estimator = "gcomp"
  )
  ## Use all 20 time points so the trapezoidal area covers [0, 20] fully.
  res <- contrast(
    fit,
    interventions = list(a1 = causatr::static(1L), a0 = causatr::static(0L)),
    times = seq_len(20L),
    type = "rmst",
    ci_method = "none"
  )
  ## Per-arm RMST at t = 20 (last grid point = the whole follow-up window).
  est <- res$estimates
  rmst_a1 <- est[get("intervention") == "a1" & get("time") == 20L, rmst_hat]
  rmst_a0 <- est[get("intervention") == "a0" & get("time") == 20L, rmst_hat]

  ## Collapse to one row per id for survRM2.
  id_dat <- dat[,
    list(
      time = if (any(get("Y") == 1L)) min(get("t")[get("Y") == 1L]) else 20L,
      status = as.integer(any(get("Y") == 1L)),
      arm = get("A")[1L]
    ),
    by = "id"
  ]

  rms <- survRM2::rmst2(
    id_dat$time,
    id_dat$status,
    id_dat$arm,
    tau = 20L
  )

  ## Extract arm-specific RMST estimates.
  km_a1 <- as.numeric(rms$RMST.arm1$rmst["Est."])
  km_a0 <- as.numeric(rms$RMST.arm0$rmst["Est."])

  ## Pooled-logistic RMST (trapezoidal) vs KM-based RMST (continuous area).
  ## Tolerance 0.05 covers the discretisation bias + finite-sample noise.
  expect_equal(rmst_a1, km_a1, tolerance = 0.05)
  expect_equal(rmst_a0, km_a0, tolerance = 0.05)
})
