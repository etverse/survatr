#' Trapezoidal cumulative RMST
#'
#' Restricted mean survival time at t* is `RMST(t*) = int_0^{t*} S(u) du`.
#' For Track A on discrete periods we approximate with the trapezoidal rule
#' on `times`, prepending `t_0 = 0` with `S(0) = 1` so the integral starts
#' at the origin:
#'
#'   RMST(t_j) = sum_{i=0}^{j-1} (t_{i+1} - t_i) * (S_i + S_{i+1}) / 2
#'
#' with `t_0 = 0`, `S_0 = 1`. This returns RMST at **every** time in
#' `times`, not just the last one -- users who want only `t* = max(times)`
#' filter after the fact.
#'
#' @param times Numeric vector, sorted ascending, first value > 0.
#' @param s_hat Numeric vector of survival estimates at each `times` entry.
#'
#' @return Numeric vector of cumulative RMSTs, same length as `times`.
#' @noRd
trapezoidal_rmst <- function(times, s_hat) {
  if (length(times) != length(s_hat)) {
    rlang::abort("`times` and `s_hat` must have equal length.")
  }
  t_full <- c(0, times)
  s_full <- c(1, s_hat)
  dt <- diff(t_full)
  ## Average of consecutive survival values: mean(S(t_i), S(t_{i+1})).
  avg_s <- (s_full[-length(s_full)] + s_full[-1L]) / 2
  cumsum(dt * avg_s)
}

#' Trapezoidal weights on a time grid
#'
#' Return the vector `w` such that the cumulative RMST at `times` equals
#' the dot product `crossprod(w_j, s_hat)` (for the final element `t*`).
#' Used by the variance engine in chunk 3 to propagate a per-individual IF
#' on `S^a_i(t)` to a per-individual IF on `RMST^a_i(t)` via `IF_RMST =
#' IF_S %*% w`.
#'
#' Not used by chunk 2's `trapezoidal_rmst()` (which just computes the
#' value); exposed now so the RMST quadrature lives in one place.
#'
#' @param times Numeric vector, sorted ascending, first value > 0.
#'
#' @return Numeric matrix of dimension `length(times) x length(times)`
#'   whose `j`th row gives the trapezoidal weights for `RMST(t_j)` on the
#'   `S(t_1), ..., S(t_{length(times)})` vector.
#' @noRd
rmst_weights <- function(times) {
  K <- length(times)
  t_full <- c(0, times)
  dt <- diff(t_full)
  W <- matrix(0, nrow = K, ncol = K)
  ## Row j holds the weights so that W[j, ] %*% s_hat = RMST(t_j).
  ## The (0, S(0) = 1) prefix contributes a constant `dt[1] / 2` that does
  ## not multiply any `s_hat` element; the variance aggregation later only
  ## needs derivatives wrt S, so the constant drops out.
  for (j in seq_len(K)) {
    for (i in seq_len(j)) {
      if (i == 1L) {
        W[j, i] <- W[j, i] + dt[1L] / 2
      } else {
        W[j, i] <- W[j, i] + dt[i - 1L] / 2
      }
      W[j, i] <- W[j, i] + dt[i] / 2
    }
  }
  W
}

#' Attach cumulative RMST to an `estimates` data.table
#'
#' When `type` is `"rmst"` or `"rmst_difference"`, replace the per-time
#' `s_hat` column with the cumulative trapezoidal integral of survival.
#' `risk_hat` is dropped (undefined for RMST). Called from
#' `contrast.survatr_fit()` before contrast assembly.
#'
#' @param estimates A `data.table` from `compute_survival_curve()` stacked
#'   across interventions.
#' @param times The user-supplied time grid.
#'
#' @return A modified `data.table` with an `rmst_hat` column replacing
#'   `s_hat` / `risk_hat`.
#' @noRd
add_rmst_to_estimates <- function(estimates, times) {
  out <- data.table::copy(estimates)
  data.table::setkeyv(out, c("intervention", "time"))
  out[,
    rmst_hat := trapezoidal_rmst(time, s_hat),
    by = "intervention"
  ]
  out[, c("s_hat", "risk_hat") := NULL]
  data.table::setcolorder(
    out,
    c("intervention", "time", "rmst_hat", "se", "ci_lower", "ci_upper", "n")
  )
  out[]
}
