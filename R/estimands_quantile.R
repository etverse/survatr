#' Solve for the q-quantile of a survival curve on the discrete grid
#'
#' @description
#' The `q`-quantile of the survival time under an intervention is the smallest
#' `t` at which the survival curve drops to (or below) `1 - q`:
#' `tau_q = inf{ t : S(t) <= 1 - q }` (the median is `q = 0.5`, i.e. the
#' crossing of `S(t) = 0.5`). Found by linear interpolation between the two
#' bracketing grid points.
#'
#' @details
#' On the discrete grid `times` with monotone non-increasing `s_hat`, let
#' `p = 1 - q`. The crossing lies in the first interval whose right endpoint is
#' at or below `p`; the left bracket is the previous grid point, or the implicit
#' origin `(0, S(0) = 1)` when the curve is already below `p` at the first grid
#' time. Linear interpolation on that interval gives `tau`.
#'
#' For the implicit-function delta-method gradient (`variance_if_survival.R`),
#' the returned list also carries the bracketing column indices `lo` / `hi`
#' (into `times`), the interpolation weight `w` (so `tau = t_lo + w (t_hi -
#' t_lo)` and the IF at `tau` is `(1 - w) IF[, lo] + w IF[, hi]`), and the local
#' curve slope `slope = (S_hi - S_lo)/(t_hi - t_lo)` (`dS/dt` at `tau`, `<= 0`).
#' `lo == 0L` encodes the origin bracket, whose IF column is the zero vector
#' (`S(0) = 1` is fixed, not estimated).
#'
#' Aborts with `survatr_quantile_unreached` when the curve never reaches `1 - q`
#' on the grid (the grid does not extend far enough to identify the quantile).
#'
#' Source: implicit-function differentiation of the survival quantile; see
#' `CHUNK_12_ESTIMANDS_QUANTILE_RMTL.md`.
#'
#' @param times Numeric grid, sorted ascending, first value `> 0`.
#' @param s_hat Numeric survival estimates at each `times` entry (monotone
#'   non-increasing in `t`).
#' @param q Scalar quantile level in `(0, 1)`.
#' @param call Caller frame for the unreached-quantile abort.
#'
#' @returns A list `list(tau, lo, hi, w, slope)`.
#' @noRd
solve_quantile <- function(times, s_hat, q, call = rlang::caller_env()) {
  p <- 1 - q
  K <- length(times)
  ## Curve never crosses below p on the grid: the q-quantile is unidentified.
  if (s_hat[K] > p) {
    rlang::abort(
      c(
        paste0(
          "The survival curve never reaches 1 - q = ",
          signif(p, 4),
          " on the time grid (min S = ",
          signif(s_hat[K], 4),
          ")."
        ),
        i = paste0(
          "The q = ",
          q,
          " quantile is beyond max(times) = ",
          times[K],
          "; extend `times` or lower `q`."
        )
      ),
      class = "survatr_quantile_unreached",
      call = call
    )
  }
  ## First grid index at or below p.
  j <- which(s_hat <= p)[1L]
  ## Left bracket: previous grid point, or the origin (0, 1) when j == 1.
  if (j == 1L) {
    t_lo <- 0
    s_lo <- 1
    lo <- 0L
  } else {
    t_lo <- times[j - 1L]
    s_lo <- s_hat[j - 1L]
    lo <- j - 1L
  }
  t_hi <- times[j]
  s_hi <- s_hat[j]
  hi <- j
  slope <- (s_hi - s_lo) / (t_hi - t_lo) ## local dS/dt at the crossing (<= 0)
  ## A flat segment (slope == 0) means S == p across the whole interval; take
  ## the right edge as the infimum and a zero gradient (SE is then degenerate;
  ## bootstrap is the documented fallback for near-flat curves).
  if (slope == 0) {
    tau <- t_hi
    w <- 1
  } else {
    ## s_lo + slope * (tau - t_lo) = p  =>  tau = t_lo + (p - s_lo)/slope.
    tau <- t_lo + (p - s_lo) / slope
    w <- (tau - t_lo) / (t_hi - t_lo)
  }
  list(tau = tau, lo = lo, hi = hi, w = w, slope = slope)
}

#' Per-individual influence function for a survival quantile
#'
#' @description
#' Propagate the `n_ids x |times|` survival IF matrix to a per-individual IF on
#' the `q`-quantile `tau_q` via the implicit-function delta method.
#'
#' @details
#' With `tau_q` solving `S(tau_q) = 1 - q`, the implicit-function theorem gives
#' `d tau_q = - dS(tau_q) / S'(tau_q)`. `dS(tau_q)` is the IF on the estimated
#' survival at the off-grid point `tau_q`, obtained by linearly interpolating
#' the two bracketing IF columns (`(1 - w) IF[, lo] + w IF[, hi]`), and
#' `S'(tau_q)` is the local grid slope (`sol$slope`). The origin bracket
#' (`lo == 0`) contributes a zero IF column (`S(0) = 1` is fixed). When the
#' slope is zero (flat curve at the crossing) the gradient is undefined; the IF
#' is returned as all-zero so the SE is `0` and the bootstrap fallback applies.
#'
#' @param if_mat The `n_ids x |times|` survival IF matrix.
#' @param sol The output of `solve_quantile()`.
#'
#' @returns Numeric vector of length `n_ids`: the per-individual IF on `tau_q`.
#' @noRd
quantile_if_vector <- function(if_mat, sol) {
  n_ids <- nrow(if_mat)
  if (sol$slope == 0) {
    return(rep(0, n_ids))
  }
  if_hi <- if_mat[, sol$hi]
  ## Origin bracket: S(0) = 1 is fixed, so its IF column is zero.
  if_lo <- if (sol$lo == 0L) rep(0, n_ids) else if_mat[, sol$lo]
  if_s_tau <- (1 - sol$w) * if_lo + sol$w * if_hi
  ## d tau = - dS(tau) / S'(tau); slope <= 0 so this flips the sign.
  -if_s_tau / sol$slope
}

#' Assemble the q-indexed quantile result (estimates + contrasts)
#'
#' @description
#' For each `(intervention, q)` solve the survival quantile from the per-
#' intervention curve; when `if_list` is supplied, also propagate the
#' implicit-function SE and Wald CI. Difference-in-quantile contrasts (e.g.
#' difference in median survival time) are built per `q` when there are at least
#' two interventions.
#'
#' @details
#' Unlike the time-indexed estimands, the quantile collapses the time axis: the
#' result is indexed by `q` (one `tau_hat` per `(intervention, q)`), so it has a
#' `q` column where the survival / RMST tables have a `time` column. SEs come
#' from `quantile_if_vector()`; the difference contrast IF is the per-individual
#' difference `IF_tau(a1) - IF_tau(ref)`.
#'
#' @param s_by_iv Named list of survival curves (numeric, ordered by `times`),
#'   one per intervention.
#' @param if_list Named list of `n_ids x |times|` IF matrices per intervention,
#'   or `NULL` for point estimates only.
#' @param times Numeric grid.
#' @param q_vec Numeric vector of quantile levels in `(0, 1)`.
#' @param iv_names Intervention names (controls result row order).
#' @param reference Reference intervention name for the contrasts, or `NULL`.
#' @param conf_level Confidence level for the Wald CI.
#' @param n_ids Number of individuals (IF normalization).
#'
#' @returns A list `list(estimates, contrasts)` of q-indexed `data.table`s.
#' @noRd
assemble_quantile_result <- function(
  s_by_iv,
  if_list,
  times,
  q_vec,
  iv_names,
  reference,
  conf_level,
  n_ids
) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  ## Stash each (iv, q) solution + per-individual IF so contrasts reuse them.
  sols <- list()
  if_tau <- list()
  est_rows <- list()
  ix <- 0L
  for (iv in iv_names) {
    s_hat <- s_by_iv[[iv]]
    for (qq in q_vec) {
      sol <- solve_quantile(times, s_hat, qq)
      key <- paste0(iv, "::", qq)
      sols[[key]] <- sol
      se <- NA_real_
      cl <- NA_real_
      cu <- NA_real_
      if (!is.null(if_list)) {
        ifv <- quantile_if_vector(if_list[[iv]], sol)
        if_tau[[key]] <- ifv
        se <- sqrt(sum(ifv^2)) / n_ids
        cl <- sol$tau - z * se
        cu <- sol$tau + z * se
      }
      ix <- ix + 1L
      est_rows[[ix]] <- data.table::data.table(
        intervention = iv,
        q = qq,
        tau_hat = sol$tau,
        se = se,
        ci_lower = cl,
        ci_upper = cu,
        n = n_ids
      )
    }
  }
  estimates <- data.table::rbindlist(est_rows)

  empty_ctr <- data.table::data.table(
    contrast = character(0),
    q = numeric(0),
    estimate = numeric(0),
    se = numeric(0),
    ci_lower = numeric(0),
    ci_upper = numeric(0)
  )
  ## A single intervention is a valid request (just the median under one arm);
  ## return the empty contrast stub rather than aborting.
  if (is.null(reference) || length(iv_names) < 2L) {
    return(list(estimates = estimates, contrasts = empty_ctr))
  }

  other <- setdiff(iv_names, reference)
  ctr_rows <- list()
  cx <- 0L
  for (a1 in other) {
    cn <- paste0(a1, " vs ", reference)
    for (qq in q_vec) {
      key_a1 <- paste0(a1, "::", qq)
      key_ref <- paste0(reference, "::", qq)
      est <- sols[[key_a1]]$tau - sols[[key_ref]]$tau
      se <- NA_real_
      cl <- NA_real_
      cu <- NA_real_
      if (!is.null(if_list)) {
        ifd <- if_tau[[key_a1]] - if_tau[[key_ref]]
        se <- sqrt(sum(ifd^2)) / n_ids
        cl <- est - z * se
        cu <- est + z * se
      }
      cx <- cx + 1L
      ctr_rows[[cx]] <- data.table::data.table(
        contrast = cn,
        q = qq,
        estimate = est,
        se = se,
        ci_lower = cl,
        ci_upper = cu
      )
    }
  }
  list(estimates = estimates, contrasts = data.table::rbindlist(ctr_rows))
}

#' Finalize a quantile contrast (point + sandwich / bootstrap variance)
#'
#' @description
#' Shared tail for the quantile estimand across all paths (gcomp / IPW / IPCW
#' main path, Track B, competing-risks all-cause). Given the per-intervention
#' survival curves and -- under sandwich -- the per-intervention IF matrices,
#' assemble the q-indexed result and fill the variance: implicit-function delta
#' for the sandwich, replicate SD for the bootstrap.
#'
#' @param fit The `survatr_fit` (needed for the bootstrap refit).
#' @param interventions Named list of interventions.
#' @param times Numeric grid.
#' @param q_vec Validated quantile levels.
#' @param reference Reference intervention name, or `NULL`.
#' @param s_by_iv Named list of survival curves (one per intervention).
#' @param if_list Named list of IF matrices (sandwich), or `NULL`.
#' @param n_ids Number of individuals.
#' @param ci_method,conf_level,n_boot,boot_ci,parallel,ncpus,seed Variance
#'   controls forwarded from `contrast()`.
#' @param call The `match.call()` for the result object.
#'
#' @returns A `survatr_result` with q-indexed `estimates` / `contrasts`.
#' @noRd
finalize_quantile <- function(
  fit,
  interventions,
  times,
  q_vec,
  reference,
  s_by_iv,
  if_list,
  n_ids,
  ci_method,
  conf_level,
  n_boot,
  boot_ci,
  parallel,
  ncpus,
  seed,
  call
) {
  res <- assemble_quantile_result(
    s_by_iv = s_by_iv,
    if_list = if_list,
    times = times,
    q_vec = q_vec,
    iv_names = names(interventions),
    reference = reference,
    conf_level = conf_level,
    n_ids = n_ids
  )
  estimates <- res$estimates
  contrasts <- res$contrasts

  if (identical(ci_method, "bootstrap")) {
    ## Resample individuals, refit, recompute the quantile per replicate; the
    ## replicate SD / percentile bands fill the SE / CI. This is the documented
    ## fallback when the curve is near-flat at tau_q and the implicit-function
    ## slope is poorly estimated.
    boot_out <- bootstrap_survival(
      fit = fit,
      interventions = interventions,
      times = times,
      type = "quantile",
      reference = reference,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed,
      q = q_vec
    )
    filled <- fill_bootstrap_ses(
      estimates = estimates,
      contrasts = contrasts,
      boot = boot_out,
      type = "quantile",
      conf_level = conf_level,
      boot_ci = boot_ci
    )
    estimates <- filled$estimates
    contrasts <- filled$contrasts
  }

  new_survatr_result(
    estimates = estimates,
    contrasts = contrasts,
    time_grid = times,
    type = "quantile",
    reference = reference,
    ci_method = ci_method,
    call = call
  )
}

#' Extract per-intervention survival curves from a time-indexed estimates table
#'
#' @description
#' Pull the `s_hat` vector (ordered by `time`) for each intervention out of the
#' chunk-2 time-indexed `estimates` table, as the input `assemble_quantile_
#' result()` needs. Track B and competing-risks all-cause curves use the same
#' `s_hat` column name, so this helper is path-agnostic.
#'
#' @param estimates A time-indexed `estimates` `data.table` with `intervention`,
#'   `time`, and `s_hat` columns.
#' @param iv_names Intervention names.
#'
#' @returns A named list of numeric survival curves, one per intervention.
#' @noRd
survival_curves_by_iv <- function(estimates, iv_names) {
  out <- lapply(iv_names, function(iv) {
    rows <- estimates[get("intervention") == iv]
    data.table::setkeyv(rows, "time")
    rows[["s_hat"]]
  })
  names(out) <- iv_names
  out
}

#' Attach cumulative RMTL to an `estimates` data.table
#'
#' @description
#' Restricted mean time **lost** up to each grid time is the complement of
#' RMST within the same window. When `type` is `"rmtl"` / `"rmtl_difference"`,
#' replace the per-time `s_hat` column with the cumulative time-lost integral.
#' `risk_hat` is dropped (undefined for RMTL). Called from
#' `contrast.survatr_fit()` before contrast assembly.
#'
#' @details
#' Restricted mean time lost at t* is
#'
#'   RMTL(t*) = t* - RMST(t*) = int_0^{t*} (1 - S(u)) du
#'
#' The `t* - RMST(t*)` form reuses `trapezoidal_rmst()` directly: the RMST
#' integral of S and the restriction time `t*` are both already on the grid,
#' so no separate quadrature of `(1 - S)` is needed. The variance reuses the
#' identical `rmst_weights()` quadratic form because the constant `t*` drops
#' out of the gradient (`d RMTL(t*) = -d RMST(t*)`), so `Var(RMTL) = Var(RMST)`.
#'
#' Source: RMTL as the complement of the restricted mean survival time;
#' see Andersen (2013) on the pseudo-observation decomposition and
#' `CHUNK_12_ESTIMANDS_QUANTILE_RMTL.md`.
#'
#' @param estimates A `data.table` from `compute_survival_curve()` stacked
#'   across interventions.
#' @param times The user-supplied time grid.
#'
#' @returns A modified `data.table` with an `rmtl_hat` column replacing
#'   `s_hat` / `risk_hat`.
#' @noRd
add_rmtl_to_estimates <- function(estimates, times) {
  out <- data.table::copy(estimates)
  data.table::setkeyv(out, c("intervention", "time"))
  ## RMTL(t) = t - RMST(t); RMST(t) is the trapezoidal area under S up to t.
  ## Computed per intervention so the cumulative integral resets at each curve.
  out[,
    rmtl_hat := time - trapezoidal_rmst(time, s_hat),
    by = "intervention"
  ]
  out[, c("s_hat", "risk_hat") := NULL]
  data.table::setcolorder(
    out,
    c("intervention", "time", "rmtl_hat", "se", "ci_lower", "ci_upper", "n")
  )
  out[]
}

#' Attach per-cause years-of-life-lost to a CIF `estimates` data.table
#'
#' @description
#' Per-cause years of life lost up to each grid time is the integral of the
#' cause-`j` cumulative incidence function. When `type = "yll"`, replace the
#' per-(intervention, cause, time) `cif_hat` column with this cumulative
#' integral. Called from `contrast_competing()` before contrast assembly.
#'
#' @details
#' Years of life lost at t* for cause `j` is
#'
#'   YLL^(j)(t*) = int_0^{t*} F^(j)(u) du
#'
#' The CIF starts at `F^(j)(0) = 0`, so the trapezoidal integral is exactly the
#' `rmst_weights()` matrix applied to the CIF vector (the `S(0) = 1` constant
#' that `trapezoidal_rmst()` carries for survival drops to `0` here). Summed over
#' causes this satisfies the identity `sum_j YLL^(j)(t*) = RMTL(t*)`, since
#' `sum_j F^(j) = 1 - S` (partition of unity) and the integral is linear.
#'
#' Source: per-cause years of life lost as the CIF integral; see Andersen (2013)
#' and `CHUNK_12_ESTIMANDS_QUANTILE_RMTL.md`.
#'
#' @param estimates A `data.table` from `compute_cif_curves()` (columns
#'   `intervention | cause | time | cif_hat | se | ci_* | n`).
#' @param times The user-supplied time grid.
#'
#' @returns A modified `data.table` with a `yll_hat` column replacing `cif_hat`.
#' @noRd
add_yll_to_estimates <- function(estimates, times) {
  out <- data.table::copy(estimates)
  data.table::setkeyv(out, c("intervention", "cause", "time"))
  ## Trapezoidal weight matrix: with F(0) = 0 the cumulative integral at each
  ## grid time is `W %*% cif_hat` (no `S(0) = 1` constant). Per (intervention,
  ## cause) so the integral resets at each CIF curve.
  weight_mat <- rmst_weights(times)
  out[,
    yll_hat := as.numeric(weight_mat %*% cif_hat),
    by = c("intervention", "cause")
  ]
  out[, "cif_hat" := NULL]
  data.table::setcolorder(
    out,
    c(
      "intervention",
      "cause",
      "time",
      "yll_hat",
      "se",
      "ci_lower",
      "ci_upper",
      "n"
    )
  )
  out[]
}
