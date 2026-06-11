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
