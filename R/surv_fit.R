#' Fit a causal survival hazard model on person-period data
#'
#' @description
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
#'   **baseline** (time-invariant) covariate adjustment set. Under
#'   `estimator = "ice"` (Track B), time-varying covariates go in
#'   `confounders_tv` instead; baseline terms here are never lagged.
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
#' @param estimator Character scalar. `"gcomp"` (pooled-logistic
#'   standardization), `"ipw"` (weighted marginal hazard MSM with stabilized
#'   density-ratio weights), or `"ice"` (Track B: longitudinal
#'   iterated-conditional-expectation hazards for a time-varying treatment).
#'   Matching is a hard reject with class `survatr_matching_rejected`
#'   pointing to `survival::coxph(..., weights = , cluster = )`.
#' @param model_fn Fitting function for the hazard model. Defaults to
#'   `stats::glm`. Accepts any function matching the `stats::glm` interface
#'   (formula, data, family, weights, ...), e.g. `mgcv::gam` with an
#'   `s(time)` term in `time_formula`.
#' @param propensity_model_fn Fitting function for the **treatment** model
#'   under `estimator = "ipw"` (ignored otherwise). Same `stats::glm`-style
#'   interface as `model_fn`. Defaults to `stats::glm`. The treatment model is
#'   fit on the baseline rows (one per id) with `confounders` as predictors;
#'   the hazard MSM then uses `A` only.
#' @param trim Numeric scalar in `(0, 1]` or `NULL` (the default). Under
#'   `estimator = "ipw"`, the per-id stabilized weights are winsorized at the
#'   `trim`-th quantile (Cole & Hernán 2008) before being broadcast onto the
#'   person-period rows. `NULL` / `1` means no truncation. The resolved fixed
#'   cutoff is reused by the sandwich variance.
#' @param confounders_tv A one-sided formula of **time-varying** confounders
#'   for Track B (`estimator = "ice"`), lag-expanded at each backward step
#'   (e.g. `~ L` builds `L + lag1_L + ...`). `NULL` (the default) means no
#'   time-varying confounders. Ignored by Track A (`gcomp` / `ipw`).
#' @param history Markov lag order for Track B. `Inf` (the default) uses the
#'   full available history (capped at `n_times - 1`); an integer restricts
#'   the lag structure (e.g. `history = 1` for first-order Markov). Ignored by
#'   Track A.
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
#' @family survatr_fit functions
#' @examples
#' # Small rectangular person-period dataset: 30 ids over 4 periods.
#' set.seed(1)
#' n_id <- 30L
#' K <- 4L
#' pp <- data.frame(
#'   id = rep(seq_len(n_id), each = K),
#'   t = rep(seq_len(K), times = n_id),
#'   A = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   Y = rbinom(n_id * K, 1L, 0.1)
#' )
#'
#' # Pooled-logistic hazard with period dummies for the baseline hazard.
#' fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
#' fit
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
  propensity_model_fn = stats::glm,
  trim = NULL,
  confounders_tv = NULL,
  history = Inf,
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
  valid_estimators <- c("gcomp", "ipw", "ice")
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
  ## Track B (ice) adds the time-varying confounders to the predictor set
  ## (Track A ignores `confounders_tv`); `time_formula` is a Track-A baseline-
  ## hazard spec and is not part of the ICE per-step formula.
  predictor_cols <- unique(c(
    outcome,
    treatment,
    id,
    time,
    all.vars(confounders),
    if (identical(estimator, "ice")) {
      all.vars(confounders_tv)
    } else {
      all.vars(time_formula)
    }
  ))
  check_no_na_in_predictors(data, predictor_cols)

  ## Outcome and censoring columns must be 0/1 indicators. `build_risk_set`
  ## and `is_uncensored()` interpret any non-binary value as if it were
  ## censored / event-positive, which produces a risk set that disagrees
  ## with the user's intent without a warning.
  check_indicator_col(data, outcome, role = "outcome", allow_na = FALSE)
  if (!is.null(censoring)) {
    check_indicator_col(data, censoring, role = "censoring", allow_na = TRUE)
  }

  check_weights(weights, nrow(data))
  check_trim(trim)

  ## Composing external (survey / design) weights with the stabilized IPW
  ## weights is a target-population transport problem (causatr handles it for
  ## scalar outcomes via a sampling block); folding it silently into the
  ## broadcast would mis-weight the MSM. It is a planned transport extension,
  ## so reject the combination upfront rather than guess the intended product.
  if (identical(estimator, "ipw") && !is.null(weights)) {
    rlang::abort(
      c(
        "External `weights` are not yet supported with `estimator = \"ipw\"`.",
        i = paste0(
          "Composing survey / external design weights with the stabilized ",
          "IPW weights (target-population transport) ships in a later chunk."
        )
      ),
      class = "survatr_ipw_external_weights"
    )
  }

  ## External / IPCW weights with Track B (ice) need the per-step weight
  ## propagation + censoring-model stacked-EE blocks of the built-in IPCW
  ## path; that ships in a later chunk. Reject upfront rather than fit an
  ## unweighted ICE chain that silently ignores the weights.
  if (identical(estimator, "ice") && !is.null(weights)) {
    rlang::abort(
      c(
        "External `weights` are not yet supported with `estimator = \"ice\"`.",
        i = paste0(
          "Weighted / IPCW longitudinal survival (per-step weight ",
          "propagation) ships in a later chunk."
        )
      ),
      class = "survatr_ice_external_weights"
    )
  }

  fit_rows <- build_risk_set(
    data = data,
    outcome = outcome,
    id = id,
    censoring = censoring
  )

  ## Estimator dispatch. gcomp fits the pooled-logistic hazard on the at-risk
  ## rows directly; ipw fits a baseline treatment model, forms stabilized
  ## weights, and fits a weighted marginal hazard MSM (alpha(t) + A); ice
  ## (Track B) fits no model here -- the per-step models are fit lazily per
  ## (intervention, horizon) inside contrast() -- and only assembles the
  ## per-step metadata (terms, families, lag order). gcomp / ipw return
  ## `model` / `family_name` / `n_fit`; ipw additionally returns the treatment
  ## models, broadcast weights, and resolved trim cutoff; ice returns
  ## `ice_details` (and `model = NULL`).
  if (identical(estimator, "ice")) {
    fit <- fit_ice_survival(
      data = data,
      fit_rows = fit_rows,
      outcome = outcome,
      treatment = treatment,
      confounders = confounders,
      confounders_tv = confounders_tv,
      id = id,
      time = time,
      history = history,
      model_fn = model_fn,
      ...
    )
    track <- "B"
  } else {
    ## Track A gcomp. A treatment that varies within id cannot be represented
    ## by the single `beta_A` of a pooled-logistic hazard MSM; warn
    ## (non-blocking) and point to the ICE estimator. Some users intentionally
    ## model a time-updated A via their own lag terms, so this is a warning, not
    ## an abort. (IPW has its own hard `survatr_ipw_time_varying_treatment`
    ## rejection downstream, so we do not also warn there.)
    if (
      identical(estimator, "gcomp") &&
        treatment_is_time_varying(data, treatment, id)
    ) {
      rlang::warn(
        c(
          paste0(
            "Treatment `",
            treatment,
            "` varies within `",
            id,
            "`, but `estimator = \"",
            estimator,
            "\"` (Track A) fits a single point-treatment hazard."
          ),
          i = "Use `estimator = \"ice\"` for a time-varying treatment."
        ),
        class = "survatr_tv_treatment_track_a"
      )
    }
    if (identical(estimator, "ipw")) {
      fit <- fit_ipw_survival(
        data = data,
        fit_rows = fit_rows,
        outcome = outcome,
        treatment = treatment,
        confounders = confounders,
        id = id,
        time = time,
        time_formula = time_formula,
        propensity_model_fn = propensity_model_fn,
        trim = trim,
        model_fn = model_fn,
        ...
      )
    } else {
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
    }
    track <- "A"
  }

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
    track = track,
    estimator = estimator,
    family = fit$family_name,
    model_fn = model_fn,
    time_formula = time_formula,
    ## For ipw, `weights` holds the per-PP-row broadcast stabilized weight;
    ## for gcomp it is the user's external weights (or NULL).
    weights = if (identical(estimator, "ipw")) fit$weights else weights,
    n_fit = fit$n_fit,
    n_total = nrow(data),
    competing = competing,
    treatment_model = fit$treatment_model,
    marginal_model = fit$marginal_model,
    trim_threshold = if (is.null(fit$trim_threshold)) {
      NA_real_
    } else {
      fit$trim_threshold
    },
    propensity_model_fn = if (identical(estimator, "ipw")) {
      propensity_model_fn
    } else {
      NULL
    },
    trim = trim,
    ## Track B (ice) metadata; NULL for Track A so the constructor defaults
    ## are unchanged.
    confounders_tv = if (identical(estimator, "ice")) confounders_tv else NULL,
    history = if (identical(estimator, "ice")) history else NULL,
    ice_details = fit$ice_details,
    call = match.call()
  )
}
