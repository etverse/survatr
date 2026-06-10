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
#' @param treatment Character scalar. Column name of the treatment. For Track A
#'   the treatment is constant within `id`. Under `estimator = "ice"` (Track B)
#'   the treatment may vary within `id` and must be **numeric** (binary, or a
#'   numeric dose entered linearly); factor / categorical (k > 2) treatments are
#'   rejected with class `survatr_ice_treatment_unsupported` (a treatment-design
#'   formula path ships in a later chunk).
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
#' @param competing Character scalar or `NULL` (the default). When non-`NULL`,
#'   activates the **competing-risks** path: `competing` names a multi-valued
#'   event-type column (`0` = no event this period, `1..J` = the cause of the
#'   event this period), and `surv_fit()` fits `J` parallel cause-specific
#'   pooled-logistic hazard models on a shared all-cause risk set. It must name
#'   the **same** column as `outcome`; competing risks are gcomp / Track A only
#'   in this release (a non-`"gcomp"` estimator, fewer than two distinct causes,
#'   or `outcome != competing` aborts with `survatr_competing_estimator` /
#'   `survatr_competing_misuse`). Cumulative-incidence functions (CIF), CIF
#'   contrasts, and all-cause survival come from `contrast()`. Fine--Gray /
#'   subdistribution hazards are out of scope (cause-specific only).
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
#'   `estimator = "ipw"`, the per-id stabilized treatment weights are
#'   winsorized at the `trim`-th quantile (Cole & Hernán 2008) before
#'   broadcast. Under IPCW (`ipcw` non-`NULL`), the same quantile is used for
#'   the **per-period** censoring weights (applied separately at each time
#'   period, targeting the heaviest late-time tails). `NULL` / `1` means no
#'   truncation. All resolved fixed cutoffs are reused by the sandwich.
#' @param confounders_tv A one-sided formula of **time-varying** confounders
#'   for Track B (`estimator = "ice"`), lag-expanded at each backward step
#'   (e.g. `~ L` builds `L + lag1_L + ...`). `NULL` (the default) means no
#'   time-varying confounders. Ignored by Track A (`gcomp` / `ipw`).
#' @param history Markov lag order for Track B. `Inf` (the default) uses the
#'   full available history (capped at `n_times - 1`); an integer restricts
#'   the lag structure (e.g. `history = 1` for first-order Markov). Ignored by
#'   Track A.
#' @param ipcw One-sided formula for the **censoring-model covariates** (e.g.
#'   `~ L1 + L2`) or `NULL` (the default, no built-in IPCW). When non-`NULL`,
#'   survatr fits a per-period censoring hazard on the person-period grid and
#'   forms stabilized **per-period cumulative** inverse-probability-of-
#'   censoring weights `W^C_{i,k}`. These are multiplied into the IPW weighted
#'   hazard MSM row weight, so the combined row weight is
#'   `w_i * W^C_{i,k}` (treatment weight × censoring weight). Requires
#'   `estimator = "ipw"` and a non-`NULL` `censoring` column; activating IPCW
#'   with any other estimator or without a censoring column aborts with a
#'   classed error. The censoring column (`censoring =`) switches from a pure
#'   row filter to the modelled path: the at-risk row set is unchanged (hazard
#'   MSM still fit on uncensored rows), but the IPCW weights reweight
#'   survivors to account for informative censoring.
#' @param censoring_model_fn Fitting function for the censoring hazard model
#'   under `ipcw` (default `stats::glm`, same `stats::glm`-style interface as
#'   `model_fn`). Ignored when `ipcw = NULL`. Stored in the fit so the
#'   bootstrap can refit the censoring model per replicate.
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
  ipcw = NULL,
  censoring_model_fn = stats::glm,
  ...
) {
  ## Reject a `MatchIt` output object passed as `data`. Users sometimes feed
  ## the `matchit` result directly; detecting it here gives a clear redirect
  ## rather than a confusing "column not found" error later.
  if (inherits(data, "matchit")) {
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

  ## Estimator gating. Matching is structurally invalid for survival (see
  ## hard-rules.md) and gets its own classed error. A `method = "matching"`
  ## mis-call (causatr-style API confusion) is also caught here. "ipw" and
  ## "ice" land on a separate classed error so downstream chunks can
  ## pattern-match when they wire them up. Anything else falls through to a
  ## generic bad-estimator error.
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
  dots <- list(...)
  if (
    identical(dots[["method"]], "matching") ||
      identical(dots[["method"]], "match")
  ) {
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

  ## Competing-risks entry point. `competing = <col>` activates the
  ## cause-specific hazards + CIF path. Two guards fire here (before touching
  ## the data) so the misuse cases keep aborting:
  ##   - CR is gcomp / Track A only this release (IPW / ICE competing risks are
  ##     later chunks): a non-gcomp estimator is a structural mismatch.
  ##   - The documented API passes the same multi-valued cause column to both
  ##     `outcome` and `competing`; a mismatch is ambiguous (which column is the
  ##     event?) so we reject rather than guess.
  is_competing <- !is.null(competing)
  if (is_competing) {
    if (!identical(estimator, "gcomp")) {
      rlang::abort(
        c(
          paste0(
            "Competing-risks (`competing =`) is supported only with ",
            "`estimator = \"gcomp\"` (Track A) in this release."
          ),
          i = paste0(
            "IPW / ICE competing-risks survival ship in later chunks. Got ",
            "`estimator = \"",
            estimator,
            "\"`."
          )
        ),
        class = "survatr_competing_estimator"
      )
    }
    if (!identical(outcome, competing)) {
      rlang::abort(
        c(
          paste0(
            "When `competing` is set it must name the same multi-valued ",
            "event-type column as `outcome` (0 = no event, 1..J = cause)."
          ),
          i = paste0(
            "Got `outcome = \"",
            outcome,
            "\"`, `competing = \"",
            competing,
            "\"`. Pass the cause column to both, or drop `competing` for the ",
            "single-event path."
          )
        ),
        class = "survatr_competing_misuse"
      )
    }
    if (!is.null(weights)) {
      rlang::abort(
        c(
          "External `weights` are not yet supported with `competing =`.",
          i = paste0(
            "Weighted / IPCW competing-risks survival ships in a later chunk."
          )
        ),
        class = "survatr_competing_weights"
      )
    }
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

  ## Outcome / censoring columns must be 0/1 indicators. `build_risk_set` and
  ## `is_uncensored()` interpret any non-binary value as if it were censored /
  ## event-positive, which produces a risk set that disagrees with the user's
  ## intent without a warning. In the competing-risks path the event column is
  ## multi-valued (0..J) by design, so it is validated by
  ## `check_competing_col()` instead and the single-event 0/1 check is skipped.
  if (is_competing) {
    causes <- check_competing_col(data, competing)
    ## Fewer than two competing causes is not a competing-risks problem -- it is
    ## the single-event path with extra ceremony. Treat it as misuse of the CR
    ## entry point rather than silently fitting one cause-specific hazard.
    if (length(causes) < 2L) {
      rlang::abort(
        c(
          paste0(
            "`competing = \"",
            competing,
            "\"` has only ",
            length(causes),
            " distinct positive cause label(s); competing risks needs >= 2."
          ),
          i = paste0(
            "With a single event type use the single-event path: drop ",
            "`competing` and pass a 0/1 `outcome`."
          )
        ),
        class = "survatr_competing_misuse"
      )
    }
  } else {
    causes <- NULL
    check_indicator_col(data, outcome, role = "outcome", allow_na = FALSE)
  }
  if (!is.null(censoring)) {
    check_indicator_col(data, censoring, role = "censoring", allow_na = TRUE)
  }

  check_weights(weights, nrow(data))
  check_trim(trim)
  check_ipcw(ipcw, estimator, censoring)

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

  ## Risk set. The single-event path keys on the 0/1 outcome; the
  ## competing-risks path keys on a derived all-cause event indicator
  ## `1{competing != 0}` so an individual leaves the risk set at the first event
  ## of *any* cause -- the shared risk set every cause-specific model fits on.
  if (is_competing) {
    data[, .survatr_any_event := as.integer(data[[competing]] != 0)]
    fit_rows <- build_risk_set(
      data = data,
      outcome = ".survatr_any_event",
      id = id,
      censoring = censoring
    )
  } else {
    fit_rows <- build_risk_set(
      data = data,
      outcome = outcome,
      id = id,
      censoring = censoring
    )
  }

  ## Estimator dispatch. Competing risks (gcomp on a multi-cause column) fits J
  ## cause-specific hazards on the shared risk set and returns
  ## `cause_models` / `causes` (with `model = NULL`). gcomp fits the
  ## pooled-logistic hazard on the at-risk rows directly; ipw fits a baseline
  ## treatment model, forms stabilized weights, and fits a weighted marginal
  ## hazard MSM (alpha(t) + A); ice (Track B) fits no model here -- the per-step
  ## models are fit lazily per (intervention, horizon) inside contrast() -- and
  ## only assembles the per-step metadata (terms, families, lag order). gcomp /
  ## ipw return `model` / `family_name` / `n_fit`; ipw additionally returns the
  ## treatment models, broadcast weights, and resolved trim cutoff; ice returns
  ## `ice_details` (and `model = NULL`).
  if (is_competing) {
    ## Competing risks is point-treatment Track A; the same single-`beta_A`
    ## caveat as gcomp applies to a within-id-varying treatment.
    if (treatment_is_time_varying(data, treatment, id)) {
      rlang::warn(
        c(
          paste0(
            "Treatment `",
            treatment,
            "` varies within `",
            id,
            "`, but competing-risks fits point-treatment cause-specific ",
            "hazards (a single `beta_A` per cause)."
          ),
          i = "Time-varying-treatment competing risks ship in a later chunk."
        ),
        class = "survatr_tv_treatment_track_a"
      )
    }
    fit <- fit_competing_risks(
      data = data,
      fit_rows = fit_rows,
      competing = competing,
      treatment = treatment,
      confounders = confounders,
      time_formula = time_formula,
      model_fn = model_fn,
      causes = causes,
      ...
    )
    track <- "A"
  } else if (identical(estimator, "ice")) {
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
      ## When IPCW is active, the censoring hazard is modelled (not just a
      ## row filter) and its cumulative inverse-probability weight multiplies
      ## the treatment weight. The IPW fit runs first (treatment weights
      ## only), then the censoring model is fit and the combined weight is
      ## passed back to the hazard MSM. Because `fit_ipw_survival()` calls
      ## `fit_hazard_gcomp()` internally, we need to intercept after the
      ## treatment weight is formed but before the MSM fit. We do this by
      ## fitting the censoring model here in surv_fit() and composing.
      if (!is.null(ipcw)) {
        fit <- fit_ipw_survival_ipcw(
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
          ipcw_formula = ipcw,
          censoring_model_fn = censoring_model_fn,
          censoring = censoring,
          ...
        )
      } else {
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
      }
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
    ## For ipw, `weights` holds the per-PP-row combined weight (treatment *
    ## IPCW when IPCW is active, or treatment only otherwise); for gcomp it
    ## is the user's external weights (or NULL).
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
    ## Competing-risks metadata; both NULL for the single-event path (the fit
    ## lists from gcomp / ipw / ice carry no `cause_models` / `causes` keys).
    cause_models = fit$cause_models,
    causes = fit$causes,
    ## IPCW metadata; populated only when `ipcw != NULL`.
    censoring_model = fit$censoring_model,
    ipcw_numerator_model = fit$ipcw_numerator_model,
    ipcw_weights = fit$ipcw_weights,
    ipcw_trim_thresholds = fit$ipcw_trim_thresholds,
    ipw_treatment_weights_pp = fit$ipw_treatment_weights_pp,
    ipcw = if (!is.null(ipcw)) ipcw else NULL,
    censoring_model_fn = if (!is.null(ipcw)) censoring_model_fn else NULL,
    call = match.call()
  )
}
