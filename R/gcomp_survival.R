#' Fit the pooled-logistic hazard model (Track A, gcomp)
#'
#' Internal fit-only helper. Assembles the model formula
#' `outcome ~ alpha(t) + A + L` from `time_formula`, `treatment`, and
#' `confounders`, subsets `data` to the at-risk rows, and calls `model_fn`
#' (default `stats::glm`) with the appropriate binomial / quasibinomial family.
#'
#' Family choice is load-bearing:
#'
#' - **Unweighted** => `stats::binomial()`. Standard score equations.
#' - **Weighted** => `stats::quasibinomial()`. Same score equations, free
#'   dispersion parameter; drops the "non-integer #successes in a binomial
#'   glm!" warning without changing the coefficient estimates.
#'
#' Swapping these families silently changes the reported SEs (dispersion
#' scaling) -- the variance engine in later chunks needs to know which family
#' was used, so the choice is recorded on the returned list as `family_name`.
#'
#' @param data Full person-period `data.table` (after risk-set construction;
#'   includes `.survatr_*` bookkeeping columns).
#' @param fit_rows Logical mask from `build_risk_set()`.
#' @param outcome,treatment Column names.
#' @param confounders RHS formula (e.g. `~ L1 + L2`) for covariate adjustment.
#' @param time_formula RHS formula for the baseline hazard `alpha(t)`.
#' @param weights External weights of length `nrow(data)`, or `NULL`.
#' @param model_fn Fitting function (e.g. `stats::glm`, `mgcv::gam`).
#' @param ... Forwarded to `model_fn`. Checked upstream to reject `na.exclude`.
#'
#' @return A list with `model` (the fitted object), `family_name` (character),
#'   and `n_fit` (integer).
#' @noRd
fit_hazard_gcomp <- function(
  data,
  fit_rows,
  outcome,
  treatment,
  confounders,
  time_formula,
  weights = NULL,
  model_fn = stats::glm,
  ...
) {
  ## Assemble RHS: baseline hazard terms first (so `alpha(t)` leads), then
  ## treatment, then confounders. Matches Hernan & Robins Ch. 17 convention.
  confounder_terms <- attr(stats::terms(confounders), "term.labels")
  time_terms <- attr(stats::terms(time_formula), "term.labels")
  rhs <- c(time_terms, treatment, confounder_terms)
  model_formula <- stats::reformulate(rhs, response = outcome)

  fit_data <- data[fit_rows]
  model_weights <- if (!is.null(weights)) weights[fit_rows] else NULL

  hazard_family <- if (!is.null(weights)) {
    stats::quasibinomial()
  } else {
    stats::binomial()
  }
  family_name <- hazard_family$family

  model <- model_fn(
    formula = model_formula,
    data = fit_data,
    family = hazard_family,
    weights = model_weights,
    ...
  )

  list(model = model, family_name = family_name, n_fit = sum(fit_rows))
}
