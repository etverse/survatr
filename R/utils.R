## Bookkeeping columns added internally during risk-set construction.
## Stripped from `fit$pp_data` before return so downstream code never sees them.
SURVATR_INTERNAL_COLS <- c(".survatr_prev_event", ".survatr_prev_cens")

#' S3 constructor for `survatr_fit`
#'
#' Internal constructor for the fit object returned by `surv_fit()`. Holds the
#' fitted hazard model, the person-period data (internal bookkeeping columns
#' stripped), and everything downstream code needs to build survival curves,
#' contrasts, diagnostics, and variance.
#'
#' @param model Fitted hazard model (e.g. a `glm` object).
#' @param pp_data Person-period `data.table` (user-facing; no `.survatr_*`
#'   internal columns).
#' @param treatment,outcome,confounders,id,time,censoring Column names / formula
#'   passed through from the user call.
#' @param time_grid Sorted unique time points present in `pp_data[[time]]`.
#' @param track Character scalar: `"A"` (point survival) or `"B"` (longitudinal
#'   ICE hazards). Currently only `"A"` ships.
#' @param estimator Character scalar: `"gcomp"`, `"ipw"`, or `"ice"`.
#' @param family Character scalar: `"binomial"` (unweighted) or
#'   `"quasibinomial"` (weighted) -- the family actually used for the GLM fit.
#' @param model_fn Fitting function the user passed (default `stats::glm`).
#' @param time_formula RHS formula for the baseline hazard `alpha(t)`.
#' @param weights External weights vector passed through, or `NULL`.
#' @param n_fit,n_total Integers: number of rows used to fit vs total PP rows.
#' @param competing Column name passed to `surv_fit()`'s `competing` argument,
#'   or `NULL`. Reserved for the cause-specific hazards + CIF path.
#' @param call The original `match.call()` of `surv_fit()`.
#'
#' @return A list of class `survatr_fit`.
#' @noRd
new_survatr_fit <- function(
  model,
  pp_data,
  treatment,
  outcome,
  confounders,
  id,
  time,
  censoring,
  time_grid,
  track,
  estimator,
  family,
  model_fn,
  time_formula,
  weights,
  n_fit,
  n_total,
  competing,
  call
) {
  structure(
    list(
      model = model,
      pp_data = pp_data,
      treatment = treatment,
      outcome = outcome,
      confounders = confounders,
      id = id,
      time = time,
      censoring = censoring,
      time_grid = time_grid,
      track = track,
      estimator = estimator,
      family = family,
      model_fn = model_fn,
      time_formula = time_formula,
      weights = weights,
      n_fit = n_fit,
      n_total = n_total,
      competing = competing,
      call = call
    ),
    class = "survatr_fit"
  )
}

#' S3 constructor for `survatr_result`
#'
#' Internal constructor for the curve-shaped result returned by
#' `contrast.survatr_fit()`. Holds the per-intervention `estimates` and
#' pairwise `contrasts` data.tables, the user time grid, and metadata
#' used by the `print` / `plot` / `tidy` methods.
#'
#' @param estimates `data.table` with columns `intervention | time |
#'   s_hat | risk_hat | ...` (or `rmst_hat` for RMST-shaped results).
#' @param contrasts `data.table` with columns `contrast | time | estimate
#'   | se | ci_lower | ci_upper`. Empty stub for curve-only `type`s.
#' @param time_grid Numeric vector (the `times` input, sorted unique).
#' @param type Contrast type.
#' @param reference Reference intervention name, or `NULL`.
#' @param ci_method `"none"` through chunk 2; `"sandwich"` / `"bootstrap"`
#'   from chunk 3 / 4.
#' @param call The original `match.call()` of `contrast.survatr_fit()`.
#'
#' @return A list of class `survatr_result`.
#' @noRd
new_survatr_result <- function(
  estimates,
  contrasts,
  time_grid,
  type,
  reference,
  ci_method,
  call
) {
  structure(
    list(
      estimates = estimates,
      contrasts = contrasts,
      time_grid = time_grid,
      type = type,
      reference = reference,
      ci_method = ci_method,
      call = call
    ),
    class = "survatr_result"
  )
}

#' Treat censoring indicator as "at-risk"
#'
#' Convention: `NA` or `0` in the censoring column means the row is at risk
#' (not censored out). Anything else (typically `1`) means the row is censored
#' at that period and must be dropped from the fit.
#'
#' Copied from causatr (`causatr:::is_uncensored`) to keep the convention
#' project-local and the error class `survatr_*`.
#'
#' @param data Person-period `data.table`.
#' @param censoring Column name of the censoring indicator, or `NULL` (no
#'   censoring column => everyone at risk).
#'
#' @return Logical vector of length `nrow(data)`.
#' @noRd
is_uncensored <- function(data, censoring) {
  if (is.null(censoring)) {
    return(rep(TRUE, nrow(data)))
  }
  cens <- data[[censoring]]
  is.na(cens) | cens == 0L
}
