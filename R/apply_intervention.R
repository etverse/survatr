#' Build the counterfactual person-period data under an intervention
#'
#' Thin wrapper around `causatr:::apply_intervention()` that (a) works on
#' person-period data (the intervention is applied row-wise to the
#' treatment column -- for baseline-constant Track A treatments this
#' broadcasts the intervened value to every period for each id, which is
#' exactly what we want), and (b) guarantees the returned `data.table` is a
#' copy, never a view into `fit$pp_data`.
#'
#' This is a deliberate use of `causatr:::` (sanctioned by
#' `.claude/hard-rules.md`: "causatr is the engine" -- we reuse its
#' intervention dispatch instead of reimplementing `static` / `shift` /
#' `scale_by` / `threshold` / `dynamic` locally).
#'
#' @param pp_data Person-period `data.table` (usually `fit$pp_data`).
#' @param treatment Character scalar -- the treatment column name.
#' @param intervention A `causatr_intervention` object from
#'   `causatr::static()` et al.
#'
#' @return A fresh `data.table` with the treatment column replaced by its
#'   intervened values.
#' @noRd
apply_intervention_pp <- function(pp_data, treatment, intervention) {
  causatr:::apply_intervention(pp_data, treatment, intervention)
}
