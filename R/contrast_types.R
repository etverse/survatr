#' Assemble the `contrasts` data.table from per-intervention estimates
#'
#' Difference and ratio contrasts (`risk_difference`, `risk_ratio`,
#' `rmst_difference`) compare each non-reference intervention against the
#' reference at every time in `estimates$time`. Curve-only types
#' (`survival`, `risk`, `rmst`) return an empty `contrasts` stub with the
#' canonical columns so downstream S3 methods see a stable shape.
#'
#' @param estimates `data.table` from `compute_survival_curve()` or
#'   `add_rmst_to_estimates()`, stacked across interventions.
#' @param type One of the six contrast types.
#' @param reference Reference intervention name, or `NULL` for curve-only
#'   types.
#' @param interventions Named list of interventions (needed for the label
#'   "<a1> vs <a0>").
#'
#' @return A `data.table` with columns
#'   `contrast | time | estimate | se | ci_lower | ci_upper`. Empty (0
#'   rows) when `type` is a curve-only type.
#' @noRd
build_contrasts <- function(estimates, type, reference, interventions) {
  empty_contrasts <- data.table::data.table(
    contrast = character(0),
    time = numeric(0),
    estimate = numeric(0),
    se = numeric(0),
    ci_lower = numeric(0),
    ci_upper = numeric(0)
  )

  if (type %in% c("survival", "risk", "rmst")) {
    return(empty_contrasts)
  }

  ref_rows <- estimates[get("intervention") == reference]
  other_names <- setdiff(names(interventions), reference)

  pieces <- lapply(other_names, function(a1_name) {
    a1_rows <- estimates[get("intervention") == a1_name]
    ## Join on time. Both tables have unique `time` entries within an
    ## intervention (setkeyv enforced upstream), so a merge is safe.
    merged <- merge(
      ref_rows,
      a1_rows,
      by = "time",
      suffixes = c(".ref", ".a1")
    )
    est <- switch(
      type,
      risk_difference = merged[["risk_hat.a1"]] - merged[["risk_hat.ref"]],
      risk_ratio = merged[["risk_hat.a1"]] / merged[["risk_hat.ref"]],
      rmst_difference = merged[["rmst_hat.a1"]] - merged[["rmst_hat.ref"]]
    )
    data.table::data.table(
      contrast = paste0(a1_name, " vs ", reference),
      time = merged$time,
      estimate = est,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  })
  data.table::rbindlist(pieces)
}
