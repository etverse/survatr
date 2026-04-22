#' Cumulative product within id, then average across ids
#'
#' Given per-row predicted hazards `h^a_{i,k}`, compute
#' `S^a_i(t) = prod_{k <= t} (1 - h^a_{i,k})` within each id and then
#' `S^a(t) = (1/n) sum_i S^a_i(t)` averaged across ids, at every user-
#' supplied time in `times`.
#'
#' Two invariants (see `.claude/hard-rules.md`):
#'
#' 1. **Cumulate within id before averaging.** The cumulative product is
#'    nonlinear; averaging hazards across ids before multiplying biases
#'    the survival estimate (Jensen's inequality).
#' 2. **Prediction covers every PP row.** The per-individual product
#'    `prod_{k <= t}` requires the counterfactual hazard at every period
#'    `k`, not only the at-risk rows from the fit step.
#'
#' @param pp_data Counterfactual person-period `data.table`.
#' @param hazards Numeric vector of per-row hazards, length `nrow(pp_data)`.
#' @param id,time Column names for the id and time variables.
#' @param times Numeric vector of user-supplied time points (sorted, unique,
#'   validated upstream).
#' @param intervention_name Character scalar -- carried through into the
#'   returned `intervention` column.
#'
#' @return A `data.table` with columns
#'   `intervention | time | s_hat | risk_hat | se | ci_lower | ci_upper | n`,
#'   one row per time in `times`. `se` / `ci_*` are `NA_real_` (filled in
#'   by later chunks).
#' @noRd
compute_survival_curve <- function(
  pp_data,
  hazards,
  id,
  time,
  times,
  intervention_name
) {
  pp <- data.table::copy(pp_data)
  pp[, .cf_hazard := hazards]
  data.table::setkeyv(pp, c(id, time))

  ## Within-id cumulative product of (1 - h). After this, `.cf_surv` at
  ## row (i, k) equals S^a_i(k).
  pp[, .cf_surv := cumprod(1 - .cf_hazard), by = c(id)]

  ## At each requested time t, average S^a_i(t) across ids that have a
  ## row at time == t. Track A data has a row for every (id, t) by
  ## construction of `prepare_pp_data()`, so this average is over all n
  ## individuals; IPW weighting enters in chunk 5.
  time_col <- time
  agg <- pp[
    get(time_col) %in% times,
    list(
      s_hat = mean(.cf_surv),
      risk_hat = 1 - mean(.cf_surv),
      n = .N
    ),
    by = c(time_col)
  ]
  data.table::setnames(agg, time_col, "time")

  agg[,
    `:=`(
      intervention = intervention_name,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  ]
  data.table::setcolorder(
    agg,
    c(
      "intervention",
      "time",
      "s_hat",
      "risk_hat",
      "se",
      "ci_lower",
      "ci_upper",
      "n"
    )
  )
  agg[]
}
