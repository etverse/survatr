## Track A g-computation sandwich vs an independent delicatessen M-estimation
## oracle. The foundational gcomp variance previously had only a closed-form
## analytic SE + empirical-SD + bootstrap check; this pins the full estimand
## surface (S^a(t), risk difference, risk ratio, RMST, RMST difference; point +
## sandwich SE) against a second implementation on a shared person-period
## fixture. Reference: data-raw/delicatessen_gcomp_survival.py (reads the same
## ipw_survival_data.csv the IPW oracle uses -- gcomp and IPW target the same
## counterfactual S^a(t) on the same confounded data).

test_that("gcomp sandwich matches the delicatessen oracle on shared data", {
  data_path <- test_path("fixtures", "python", "ipw_survival_data.csv")
  ref_path <- test_path("fixtures", "python", "gcomp_survival_delicatessen.csv")
  skip_if(!file.exists(ref_path))

  df <- data.table::fread(data_path)
  deli <- data.table::fread(ref_path)
  fit <- surv_fit(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    time_formula = ~ factor(t)
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  tms <- c(2L, 3L, 4L, 5L)

  sv <- contrast(
    fit,
    ivs,
    times = tms,
    type = "survival",
    ci_method = "sandwich"
  )
  rd <- contrast(
    fit,
    ivs,
    times = tms,
    type = "risk_difference",
    reference = "a1",
    ci_method = "sandwich"
  )
  rr <- contrast(
    fit,
    ivs,
    times = tms,
    type = "risk_ratio",
    reference = "a1",
    ci_method = "sandwich"
  )
  rm <- contrast(fit, ivs, times = tms, type = "rmst", ci_method = "sandwich")
  rmd <- contrast(
    fit,
    ivs,
    times = tms,
    type = "rmst_difference",
    reference = "a1",
    ci_method = "sandwich"
  )

  ## Map each delicatessen quantity to its survatr (estimate, se).
  get_one <- function(q, t) {
    if (q == "S1" || q == "S0") {
      a <- if (q == "S1") "a1" else "a0"
      r <- sv$estimates[get("intervention") == a & get("time") == t]
      return(c(r$s_hat, r$se))
    }
    if (q == "RMST1" || q == "RMST0") {
      a <- if (q == "RMST1") "a1" else "a0"
      r <- rm$estimates[get("intervention") == a & get("time") == t]
      return(c(r$rmst_hat, r$se))
    }
    cn <- "a0 vs a1"
    tbl <- switch(q, RD = rd, RR = rr, RMSTdiff = rmd)
    r <- tbl$contrasts[get("contrast") == cn & get("time") == t]
    c(r$estimate, r$se)
  }

  for (i in seq_len(nrow(deli))) {
    q <- deli$quantity[i]
    sv_vals <- get_one(q, deli$time[i])
    expect_equal(sv_vals[1L], deli$estimate[i], tolerance = 1e-3)
    ## risk-ratio SE is reported natural-scale (RR * se(log RR)) in survatr vs
    ## a direct delta on RR in the oracle; identical to first order, so a hair
    ## looser. All other SEs match the M-estimation sandwich tightly.
    se_tol <- if (q == "RR") 1e-2 else 1e-3
    expect_equal(sv_vals[2L], deli$se[i], tolerance = se_tol)
  }
})

test_that("rmtl sandwich matches the delicatessen oracle (RMTL = t - RMST)", {
  ## RMTL has no separate delicatessen estimating equation -- it IS the RMST
  ## complement, RMTL(t) = t - RMST(t) with Var(RMTL) = Var(RMST). So the
  ## delicatessen RMST point/SE in the shared fixture pin RMTL directly: the
  ## implied oracle RMTL is `time - RMST_deli` with the identical SE.
  data_path <- test_path("fixtures", "python", "ipw_survival_data.csv")
  ref_path <- test_path("fixtures", "python", "gcomp_survival_delicatessen.csv")
  skip_if(!file.exists(ref_path))

  df <- data.table::fread(data_path)
  deli <- data.table::fread(ref_path)
  fit <- surv_fit(
    df,
    outcome = "Y",
    treatment = "A",
    confounders = ~L,
    id = "id",
    time = "t",
    time_formula = ~ factor(t)
  )
  ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
  tms <- c(2L, 3L, 4L, 5L)
  rl <- contrast(fit, ivs, times = tms, type = "rmtl", ci_method = "sandwich")

  ## Each arm's RMTL vs the implied oracle (t - RMST_deli, with se_deli).
  rmst_rows <- deli[deli$quantity %in% c("RMST1", "RMST0"), ]
  for (i in seq_len(nrow(rmst_rows))) {
    a <- if (rmst_rows$quantity[i] == "RMST1") "a1" else "a0"
    tt <- rmst_rows$time[i]
    sv <- rl$estimates[get("intervention") == a & get("time") == tt]
    expect_equal(sv$rmtl_hat, tt - rmst_rows$estimate[i], tolerance = 1e-3)
    expect_equal(sv$se, rmst_rows$se[i], tolerance = 1e-3)
  }
})
