#' Predict per-row hazard on counterfactual person-period data
#'
#' Calls `stats::predict(model, newdata = pp_cf_data, type = "response")` to
#' get the predicted hazard at every row. Unlike the fit step, prediction
#' runs over **all** rows -- the cumulative product `S^a_i(t) = prod_{k <= t}
#' (1 - h^a_{i,k})` needs `h^a_{i,k}` at every period `k <= t`, whether or
#' not individual `i` was still at risk in reality. Restricting prediction
#' to `fit_rows` would truncate the cumulative product at the observed
#' event / censor time and bias the estimator downward (missing the
#' counterfactual tail).
#'
#' @param model A fitted hazard model (typically `fit$model`, a `glm` or
#'   compatible object).
#' @param pp_cf_data Counterfactual person-period `data.table` from
#'   `apply_intervention_pp()`.
#'
#' @return Numeric vector of predicted hazards, length `nrow(pp_cf_data)`.
#' @noRd
predict_hazard_pp <- function(model, pp_cf_data) {
  unname(stats::predict(model, newdata = pp_cf_data, type = "response"))
}
