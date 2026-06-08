## Bookkeeping columns added internally during risk-set construction.
## Stripped from `fit$pp_data` before return so downstream code never sees them.
## `.survatr_any_event` / `.survatr_cause_event` are the transient all-cause and
## per-cause 0/1 event indicators derived from the multi-valued `competing`
## column in the cause-specific hazards path; like the lagged cumsums they are
## stripped before return.
SURVATR_INTERNAL_COLS <- c(
  ".survatr_prev_event",
  ".survatr_prev_cens",
  ".survatr_any_event",
  ".survatr_cause_event"
)

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
#'   or `NULL`. When non-`NULL` the fit is a cause-specific hazards + CIF
#'   competing-risks fit: `model` is `NULL` and the per-cause hazard models live
#'   in `cause_models`.
#' @param cause_models Named list of fitted cause-specific hazard models (one
#'   per competing event type, named by the integer cause label), or `NULL`
#'   (single-event fit). Each is a pooled-logistic hazard model on the shared
#'   all-cause risk set with the cause-`j` event indicator as response.
#' @param causes Integer vector of the competing event-type labels (`1..J`)
#'   present in the data, or `NULL` (single-event fit). Names of `cause_models`
#'   are `as.character(causes)`.
#' @param treatment_model The full `A ~ L` `causatr_treatment_model` under
#'   `estimator = "ipw"`, or `NULL` (gcomp). Used by the IPW sandwich and
#'   `diagnose()`.
#' @param marginal_model The marginal `A ~ 1` numerator
#'   `causatr_treatment_model` under `estimator = "ipw"`, or `NULL`.
#' @param trim_threshold Numeric scalar: the fixed weight-winsorization cutoff
#'   used in the IPW fit (so the sandwich clips at the same value), or
#'   `NA_real_` when no truncation / not IPW.
#' @param propensity_model_fn Fitting function used for the treatment model
#'   under `estimator = "ipw"`, or `NULL`. Stored so the bootstrap can refit the
#'   treatment model identically per replicate.
#' @param trim The `trim` quantile level passed to `surv_fit()` (not the
#'   resolved cutoff), or `NULL`. Stored so the bootstrap re-estimates the
#'   winsorization per replicate.
#' @param confounders_tv Time-varying confounders formula under
#'   `estimator = "ice"` (Track B), or `NULL` (Track A). Threaded to the
#'   bootstrap refit and the contrast ICE path.
#' @param history Markov lag order under `estimator = "ice"`, or `NULL`
#'   (Track A).
#' @param ice_details Per-step ICE metadata (`build_ice_surv_details()` output)
#'   under `estimator = "ice"`, or `NULL`. Lets `contrast()` / the bootstrap
#'   reconstruct the ICE machinery without re-parsing.
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
  call,
  treatment_model = NULL,
  marginal_model = NULL,
  trim_threshold = NA_real_,
  propensity_model_fn = NULL,
  trim = NULL,
  confounders_tv = NULL,
  history = NULL,
  ice_details = NULL,
  cause_models = NULL,
  causes = NULL
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
      treatment_model = treatment_model,
      marginal_model = marginal_model,
      trim_threshold = trim_threshold,
      propensity_model_fn = propensity_model_fn,
      trim = trim,
      confounders_tv = confounders_tv,
      history = history,
      ice_details = ice_details,
      cause_models = cause_models,
      causes = causes,
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
