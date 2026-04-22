#' Fit a causal survival hazard model on person-period data
#'
#' Fit-only entry point for survatr. Builds the risk set and fits the
#' pooled-logistic discrete-time hazard model
#' `logit h(t | A, L) = alpha(t) + beta_A A + beta_L L` on the at-risk
#' person-period rows. Survival curves, risk / RMST contrasts, and variance
#' live in `contrast()` (time-indexed curve-shaped result). This two-step
#' split lets the user fit the hazard model once and cheaply contrast many
#' interventions on top.
#'
#' @param data A person-period (long) `data.frame` or `data.table`,
#'   **rectangular** across ids: every unique `id` must have one row at
#'   every unique `time` value. Ragged PP (ids dropped post-event /
#'   post-censor) is rejected with class `survatr_ragged_pp`. Pad ragged
#'   data before calling `surv_fit()` by appending rows with
#'   `outcome = 0` and (if used) `censoring = 1` so the risk-set builder
#'   drops them from the fit. Wide data (one row per id across a
#'   multi-period study) must be reshaped with
#'   `causatr::to_person_period()`.
#' @param outcome Character scalar. Column name of the event indicator
#'   (`1` = event at this period, `0` = no event). Must be in `data`.
#' @param treatment Character scalar. Column name of the (baseline, point)
#'   treatment. For Track A the treatment is constant within `id`.
#' @param confounders A one-sided formula (e.g. `~ L1 + L2`) describing the
#'   covariate adjustment set.
#' @param id Character scalar. Column name of the individual identifier.
#' @param time Character scalar. Column name of the discrete period index
#'   (integer-valued; sorted within `id`).
#' @param censoring Character scalar or `NULL`. Column name of the censoring
#'   indicator (`1` = censored at this period). `NA` or `0` means uncensored.
#'   When `NULL`, every uncensored period is treated as at-risk until the
#'   first event.
#' @param competing Reserved for the cause-specific hazards + CIF path.
#'   Passing anything other than `NULL` is an error in the current release --
#'   competing risks have a dedicated entry point rather than being silently
#'   folded into `surv_fit()` (which would fit a biased cause-deleted
#'   hazard).
#' @param time_formula One-sided formula for the baseline hazard
#'   `alpha(t)`. Defaults to `~ splines::ns(time, 4)` (4 df natural spline
#'   on the time variable). Pass `~ factor(time)` for period dummies or
#'   `~ 1` for a time-constant hazard.
#' @param weights Optional numeric vector of external weights, length
#'   `nrow(data)`. When supplied, the hazard model is fit with
#'   `stats::quasibinomial()` rather than `stats::binomial()` (same score
#'   equations, free dispersion -- drops the "non-integer #successes"
#'   warning). The variance engine in later chunks reads the family from
#'   `fit$family` to pick the right dispersion.
#' @param estimator Character scalar. Currently `"gcomp"` (pooled-logistic
#'   standardization). `"ipw"` and `"ice"` ship in later tracks and are
#'   rejected here. Matching is a hard reject with class
#'   `survatr_matching_rejected` pointing to
#'   `survival::coxph(..., weights = , cluster = )`.
#' @param model_fn Fitting function for the hazard model. Defaults to
#'   `stats::glm`. Accepts any function matching the `stats::glm` interface
#'   (formula, data, family, weights, ...), e.g. `mgcv::gam` with an
#'   `s(time)` term in `time_formula`.
#' @param ... Forwarded to `model_fn`. `na.action = na.exclude` is rejected
#'   with class `survatr_bad_na_action` -- `na.exclude` pads working
#'   residuals with `NA`s while `model.matrix()` drops them, which silently
#'   corrupts the sandwich variance downstream.
#'
#' @return An object of class `survatr_fit` holding the fitted hazard model,
#'   the person-period data (internal `.survatr_*` columns stripped), the
#'   time grid, and metadata needed by `contrast()` and `diagnose()`.
#'
#' @seealso `causatr::to_person_period()` for reshaping wide data.
#'
#' @export
surv_fit <- function(
  data,
  outcome,
  treatment,
  confounders,
  id,
  time,
  censoring = NULL,
  competing = NULL,
  time_formula = ~ splines::ns(time, 4),
  weights = NULL,
  estimator = "gcomp",
  model_fn = stats::glm,
  ...
) {
  ## Estimator gating. Matching is structurally invalid for survival (see
  ## hard-rules.md) and gets its own classed error. "ipw" and "ice" land on
  ## a separate classed error so downstream chunks can pattern-match when
  ## they wire them up. Anything else falls through to a generic bad-estimator
  ## error.
  if (identical(estimator, "matching") || identical(estimator, "match")) {
    rlang::abort(
      c(
        "Matching + survival is out of scope for survatr.",
        i = paste0(
          "Use `survival::coxph(..., weights = match_weights, ",
          "cluster = subclass)` directly on the `MatchIt` output."
        )
      ),
      class = "survatr_matching_rejected"
    )
  }
  valid_estimators <- "gcomp"
  if (!isTRUE(estimator %in% valid_estimators)) {
    rlang::abort(
      paste0(
        "`estimator` must be one of: ",
        paste0("\"", valid_estimators, "\"", collapse = ", "),
        ". Got \"",
        estimator,
        "\"."
      ),
      class = "survatr_bad_estimator"
    )
  }

  check_dots_na_action(...)

  if (!is.null(competing)) {
    rlang::abort(
      c(
        paste0(
          "Competing-risks survival is not handled by `surv_fit()` -- the ",
          "`competing` argument is reserved for a dedicated cause-specific ",
          "hazards + cumulative-incidence path."
        ),
        i = paste0(
          "Passing `competing = NULL` (the default) to `surv_fit()` with ",
          "competing events present in `data[[",
          outcome,
          "]]` would fit a ",
          "biased cause-deleted hazard."
        )
      ),
      class = "survatr_competing_misuse"
    )
  }

  check_reserved_cols(data)

  data <- prepare_pp_data(
    data = data,
    outcome = outcome,
    treatment = treatment,
    id = id,
    time = time,
    censoring = censoring
  )

  ## Reject NA in any column that feeds into the hazard-model formula
  ## or the counterfactual prediction. `censoring` is excluded because
  ## NA there carries "uncensored" semantics via `is_uncensored()`.
  predictor_cols <- unique(c(
    outcome,
    treatment,
    id,
    time,
    all.vars(confounders),
    all.vars(time_formula)
  ))
  check_no_na_in_predictors(data, predictor_cols)

  check_weights(weights, nrow(data))

  fit_rows <- build_risk_set(
    data = data,
    outcome = outcome,
    id = id,
    censoring = censoring
  )

  fit <- fit_hazard_gcomp(
    data = data,
    fit_rows = fit_rows,
    outcome = outcome,
    treatment = treatment,
    confounders = confounders,
    time_formula = time_formula,
    weights = weights,
    model_fn = model_fn,
    ...
  )

  ## Strip internal bookkeeping before returning. Downstream code
  ## (`contrast()`, `diagnose()`) rebuilds the risk set from scratch from
  ## the user-facing columns; it must not see `.survatr_*` in `fit$pp_data`.
  internal <- intersect(SURVATR_INTERNAL_COLS, names(data))
  if (length(internal) > 0L) {
    data[, (internal) := NULL]
  }

  new_survatr_fit(
    model = fit$model,
    pp_data = data,
    treatment = treatment,
    outcome = outcome,
    confounders = confounders,
    id = id,
    time = time,
    censoring = censoring,
    time_grid = sort(unique(data[[time]])),
    track = "A",
    estimator = "gcomp",
    family = fit$family_name,
    model_fn = model_fn,
    time_formula = time_formula,
    weights = weights,
    n_fit = fit$n_fit,
    n_total = nrow(data),
    competing = competing,
    call = match.call()
  )
}
