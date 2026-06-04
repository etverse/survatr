#' Assemble the Track B (ICE survival) fit metadata
#'
#' @description
#' Fit-only entry for `estimator = "ice"`, called by `surv_fit()`. No model is
#' fit here -- the per-step ICE models are fit lazily per (intervention,
#' horizon) inside `contrast()`. This helper only validates the treatment
#' structure and assembles the per-step metadata (`build_ice_surv_details()`)
#' that `contrast()` and the bootstrap reuse.
#'
#' A treatment that is constant within id is **valid** for ICE (the strategy
#' is applied to a static treatment); we `inform()` once that Track A
#' (`estimator = "gcomp"`) is cheaper for that case, but do not abort or
#' silently switch tracks -- they are different estimators.
#'
#' @param data Person-period `data.table` (after risk-set construction).
#' @param fit_rows Logical at-risk mask from `build_risk_set()`.
#' @param outcome,treatment Column names.
#' @param confounders Baseline confounder formula.
#' @param confounders_tv Time-varying confounder formula, or `NULL`.
#' @param id,time Column names.
#' @param history Markov lag order (`Inf` = full history).
#' @param model_fn Per-step fitting function.
#' @param ... Captured and replayed per step via `ice_fit_step()`.
#'
#' @return A list with `model` (`NULL`), `family_name` (`"binomial"`, the
#'   horizon-step family), `n_fit` (at-risk row count), and `ice_details`.
#' @noRd
fit_ice_survival <- function(
  data,
  fit_rows,
  outcome,
  treatment,
  confounders,
  confounders_tv,
  id,
  time,
  history,
  model_fn = stats::glm,
  ...
) {
  ## Track B models the treatment as a numeric main effect (binary is the
  ## special case; a numeric dose enters linearly, the user's modelling
  ## choice). A factor / character treatment is NOT supported here: the
  ## intervention sets a numeric counterfactual value, which collides with a
  ## factor column (a cryptic data.table error), and nominal categories need a
  ## `treatment_form = ~ factor(A)` design that Track B does not yet thread.
  ## Reject with a classed error pointing at the future extended-types chunk.
  ## Issue #2, 2026-06-03 critical review (/tmp/survatr_repro_cat.R).
  if (!is.numeric(data[[treatment]])) {
    rlang::abort(
      c(
        paste0(
          "`estimator = \"ice\"` requires a numeric `treatment` (binary or a ",
          "numeric dose modelled linearly); got a ",
          class(data[[treatment]])[1L],
          "."
        ),
        i = paste0(
          "Categorical (k > 2) treatments via `factor()` / a treatment design ",
          "formula ship in a later chunk; recode to a numeric column for now."
        )
      ),
      class = "survatr_ice_treatment_unsupported"
    )
  }

  if (!treatment_is_time_varying(data, treatment, id)) {
    rlang::inform(
      c(
        paste0(
          "`estimator = \"ice\"` with a treatment constant within `",
          id,
          "` is valid (Track B applies the strategy to a static treatment)."
        ),
        i = paste0(
          "For a point (baseline-constant) treatment, ",
          "`estimator = \"gcomp\"` (Track A) is cheaper."
        )
      ),
      class = "survatr_ice_static_treatment",
      .frequency = "once",
      .frequency_id = "survatr_ice_static_treatment"
    )
  }

  ice_details <- build_ice_surv_details(
    confounders = confounders,
    confounders_tv = confounders_tv,
    treatment = treatment,
    time_grid = sort(unique(data[[time]])),
    history = history,
    model_fn = model_fn,
    dots = list(...),
    weights = NULL
  )

  list(
    model = NULL,
    family_name = ice_details$family_outcome$family,
    n_fit = sum(fit_rows),
    ice_details = ice_details
  )
}

#' Is the treatment time-varying within individuals?
#'
#' @description
#' Detect whether the treatment column takes more than one distinct value
#' within any individual -- the signature of a time-varying treatment that
#' Track A's single-`beta_A` pooled-logistic hazard cannot represent and that
#' Track B (ICE) is built for.
#'
#' @param data Person-period `data.table`.
#' @param treatment Treatment column name.
#' @param id Individual identifier column name.
#'
#' @return `TRUE` if at least one id has >1 distinct treatment value.
#' @noRd
treatment_is_time_varying <- function(data, treatment, id) {
  per_id <- data[,
    list(.nu = data.table::uniqueN(get(treatment))),
    by = c(id)
  ]
  any(per_id$.nu > 1L)
}

#' Track B contrast driver (longitudinal ICE-hazard survival)
#'
#' @description
#' The `contrast.survatr_fit()` branch for `fit$track == "B"`. Builds the
#' per-intervention survival curves via the ICE backward sweep
#' (`compute_ice_survival_curve()`), then reuses the Track A estimand
#' machinery wholesale: `build_contrasts()` for pairwise contrasts,
#' `add_rmst_to_estimates()` for RMST, `fill_sandwich_ses()` for the
#' delta-method CIs (fed the ICE influence-function matrices), and
#' `bootstrap_survival()` / `fill_bootstrap_ses()` for the empirical bootstrap
#' (which refits the ICE chain per replicate).
#'
#' @inheritParams contrast.survatr_fit
#' @param call The `match.call()` captured in `contrast.survatr_fit()` for the
#'   result's `call` slot.
#'
#' @return A `survatr_result`.
#' @noRd
contrast_track_b <- function(
  fit,
  interventions,
  times,
  type,
  reference,
  ci_method,
  conf_level,
  n_boot,
  boot_ci,
  parallel,
  ncpus,
  seed,
  call
) {
  ## Track B v1 supports only deterministic plug-in interventions (`static`,
  ## `shift`, `scale_by`, `threshold`, `dynamic`). Two families are rejected
  ## upfront:
  ##   - `ipsi()` reweights the propensity (weight-path estimand, IPW-only);
  ##     there is no counterfactual treatment to plug into the ICE recursion.
  ##   - Stochastic interventions need Monte-Carlo averaging at the cumulative-
  ##     product level (Jensen-safe) and an MC sensitivity chain.
  ## Both ship in later chunks; reject here rather than (for stochastic)
  ## silently return a single random draw or (for ipsi) surface causatr's
  ## generic g-computation abort.
  iv_types <- vapply(
    interventions,
    function(iv) if (is.null(iv$type)) NA_character_ else iv$type,
    character(1L)
  )
  is_stoch <- vapply(
    interventions,
    causatr:::has_stochastic_component,
    logical(1L)
  )
  if (any(iv_types == "ipsi", na.rm = TRUE) || any(is_stoch)) {
    rlang::abort(
      c(
        paste0(
          "Only deterministic interventions (`static`, `shift`, `scale_by`, ",
          "`threshold`, `dynamic`) are supported for Track B (ice)."
        ),
        i = paste0(
          "Incremental-propensity (`ipsi()`) and stochastic / Monte-Carlo ",
          "survival interventions ship in later chunks."
        )
      ),
      class = "survatr_ice_intervention_deferred"
    )
  }

  details <- fit$ice_details
  ## Lags + risk set, once (intervention-agnostic).
  base <- prepare_track_b_base(fit, details)
  n_ids <- length(unique(as.character(base$data_lag[[fit$id]])))

  ## Per-intervention survival curves (+ stashed per-horizon ice_results for
  ## the sandwich bridge).
  per_iv <- lapply(names(interventions), function(nm) {
    compute_ice_survival_curve(base, details, interventions[[nm]], times, nm)
  })
  names(per_iv) <- names(interventions)

  estimates <- data.table::rbindlist(lapply(per_iv, function(x) x$estimates))
  if (type %in% c("rmst", "rmst_difference")) {
    estimates <- add_rmst_to_estimates(estimates, times)
  }
  contrasts <- build_contrasts(estimates, type, reference, interventions)

  if (identical(ci_method, "sandwich")) {
    ## Stacked-EE ICE sandwich via causatr's chain, reused per horizon. The
    ## minimal causatr_fit + per-horizon ice_results feed
    ## `causatr:::variance_if_ice_one()`; the resulting n x |times| IF matrix
    ## drops straight into the shared `fill_sandwich_ses()`.
    min_fit <- build_min_causatr_fit_b(fit, details, base$data_lag)
    ## Target = the at-risk-at-baseline standardisation population (entry-
    ## censored ids carry NA `pseudo_final` and must be excluded from Channel-1,
    ## else `mu_hat` is NA). Aligned to the first-period id order `ice_if_setup`
    ## uses. `n_ids` stays the FULL count: `ice_if_setup` scales Channel 1 by
    ## `n / n_target`, so `crossprod(IF) / n_ids^2` is the correct mean variance.
    ## Issue #1, 2026-06-03 critical review.
    first_t <- details$time_points[1]
    first_mask <- base$data_lag[[fit$time]] == first_t
    target <- base$fit_rows[first_mask]
    ## Per-period event indicators feed the chain's `(1 - D_k)` survival
    ## failure carry-forward factor.
    event_by_step <- build_event_by_step(
      data = base$data_lag,
      time_points = details$time_points,
      id_col = fit$id,
      time_col = fit$time,
      outcome = fit$outcome
    )
    if_list <- lapply(names(interventions), function(nm) {
      compute_ice_survival_if_matrix(
        min_fit = min_fit,
        ice_results = per_iv[[nm]]$ice_results,
        times = times,
        target = target,
        event_by_step = event_by_step
      )
    })
    names(if_list) <- names(interventions)
    filled <- fill_sandwich_ses(
      estimates = estimates,
      contrasts = contrasts,
      if_list = if_list,
      type = type,
      reference = reference,
      times = times,
      conf_level = conf_level,
      n_ids = n_ids
    )
    estimates <- filled$estimates
    contrasts <- filled$contrasts
  } else if (identical(ci_method, "bootstrap")) {
    ## Resample ids, refit the ICE chain per replicate (lags rebuilt inside
    ## surv_fit), recompute curves / contrasts. `bootstrap_survival()` is
    ## track-agnostic: it refits via `surv_fit(estimator = fit$estimator)`,
    ## which dispatches to Track B and threads `confounders_tv` / `history`.
    boot_out <- bootstrap_survival(
      fit = fit,
      interventions = interventions,
      times = times,
      type = type,
      reference = reference,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed
    )
    filled <- fill_bootstrap_ses(
      estimates = estimates,
      contrasts = contrasts,
      boot = boot_out,
      type = type,
      conf_level = conf_level,
      boot_ci = boot_ci
    )
    estimates <- filled$estimates
    contrasts <- filled$contrasts
  }

  new_survatr_result(
    estimates = estimates,
    contrasts = contrasts,
    time_grid = times,
    type = type,
    reference = reference,
    ci_method = ci_method,
    call = call
  )
}

#' Survival-tail pseudo-outcome for the ICE backward step
#'
#' @description
#' Construct the per-row survival-tail pseudo-outcome used as the response of
#' the quasibinomial backward-step model in Track B. For an individual at risk
#' at period `k` with all-cause event indicator `D_k` and next-step predicted
#' tail risk `q_next` (the counterfactual `1 - Ŝ^d_{k+1:τ}` carried back from
#' the later step), the pseudo-outcome is
#' \deqn{\tilde Y_k = D_k + (1 - D_k)\,q_{next} = 1 - (1 - D_k)(1 - q_{next}).}
#'
#' @details
#' This is the one piece causatr's plain ICE engine does **not** supply:
#' `causatr:::ice_iterate()` regresses the raw next-step prediction `q_next`
#' directly (single terminal outcome, no failure carry-forward). For
#' time-to-event data the failure mass at period `k` must be carried forward,
#' which is exactly the `D_k` term. The `ifelse()` form is used rather than the
#' algebraically identical `D_k + (1 - D_k) * q_next` so a failed-at-`k`
#' individual whose `q_next` is `NA` (they are not at risk at `k+1`, so they
#' never received a next-step prediction) yields `1`, not `0 * NA = NaN`.
#'
#' @param d_k Numeric/integer vector of all-cause event indicators (0/1) at the
#'   current period, one per at-risk-at-`k` row.
#' @param q_next Numeric vector of next-step predicted tail risks aligned to
#'   `d_k`; may contain `NA` for rows that failed at `k`.
#'
#' @return Numeric vector of survival-tail pseudo-outcomes in `[0, 1]` (with
#'   `NA` retained only for survived-but-unpredictable rows -- e.g. censored at
#'   `k+1` -- which the caller drops from the fit).
#' @noRd
ice_surv_pseudo <- function(d_k, q_next) {
  ## Failed at k -> deterministic 1 (NA-safe); survived -> the carried tail.
  ifelse(d_k == 1, 1, q_next)
}

#' Build per-step metadata for Track B ICE survival
#'
#' @description
#' Resolve the survival-aware mirror of `causatr::fit_ice()`'s `details`: the
#' time grid, baseline / time-varying confounder term labels, the lag (Markov)
#' order, effect-modification info, treatment design terms, the per-step
#' families, and the user's `model_fn` / `...`. Computed from the fit metadata
#' (no data scan beyond the time grid) so it can be stashed on the
#' `survatr_fit` and reused by `contrast()` and the bootstrap.
#'
#' The per-step **link forcing** lives here and is load-bearing
#' (`.claude/hard-rules.md`): `family_outcome = binomial()` for the horizon
#' (final) step where the response is the 0/1 event indicator, and
#' `family_pseudo = quasibinomial()` for earlier steps where the survival-tail
#' pseudo-outcome is a fitted probability in `[0, 1]` (the fractional-logit
#' approach of Papke & Wooldridge 1996, applied to ICE by Zivich et al. 2024).
#'
#' @param confounders Baseline (time-invariant) confounder formula.
#' @param confounders_tv Time-varying confounder formula, or `NULL`.
#' @param treatment Treatment column name.
#' @param time_grid Sorted unique time points (`fit$time_grid`).
#' @param history Markov lag order (`Inf` = full history, capped at
#'   `n_times - 1`).
#' @param model_fn Per-step fitting function (e.g. `stats::glm`).
#' @param dots User `...` captured at `surv_fit()` time, replayed per step.
#' @param weights External weight vector or `NULL` (Track B v1: always `NULL`).
#'
#' @return A named list of ICE survival details (see body) with the same field
#'   names `causatr:::ice_iterate()` / `variance_if_ice_one()` read off a
#'   `causatr_fit$details` plus the survival per-step families.
#' @noRd
build_ice_surv_details <- function(
  confounders,
  confounders_tv,
  treatment,
  time_grid,
  history,
  model_fn,
  dots,
  weights = NULL
) {
  time_points <- sort(unique(time_grid))
  n_times <- length(time_points)

  ## Baseline confounders enter as term labels (preserving `I(age^2)`,
  ## `factor(x)`); they are time-invariant and never lagged.
  baseline_terms <- attr(stats::terms(confounders), "term.labels")

  ## Effect-modification (`A:modifier`) terms are auto lag-expanded by
  ## `ice_build_formula()`; parse them from the *baseline* confounders the
  ## same way causatr does (treatment interactions live there).
  em_info <- causatr:::parse_effect_mod(confounders, treatment)

  ## Time-varying confounders in two forms: plain names (for lag-column
  ## creation + existence checks) and term labels (preserving transforms
  ## like `ns(L, 3)` so lag expansion yields `ns(lag1_L, 3)`).
  tv_vars <- if (!is.null(confounders_tv)) {
    all.vars(confounders_tv)
  } else {
    character(0)
  }
  tv_terms <- if (!is.null(confounders_tv)) {
    attr(stats::terms(confounders_tv), "term.labels")
  } else {
    character(0)
  }

  ## Bare-numeric treatment main effect (survatr Track B does not expose a
  ## `treatment_form`; the intervention always sets the numeric column).
  treatment_terms <- treatment

  ## Resolve the Markov order: `Inf` -> full history (`n_times - 1`, no lag at
  ## the first period); an integer caps the lag structure.
  max_lag <- if (is.infinite(history)) {
    n_times - 1L
  } else {
    as.integer(history)
  }

  list(
    time_points = time_points,
    n_times = n_times,
    baseline_terms = baseline_terms,
    tv_vars = tv_vars,
    tv_terms = tv_terms,
    max_lag = max_lag,
    em_info = em_info,
    treatment_terms = treatment_terms,
    model_fn = model_fn,
    ## Per-step link forcing (load-bearing -- see @description).
    family_outcome = stats::binomial(),
    family_pseudo = stats::quasibinomial(),
    dots = dots,
    weights = weights,
    stratified = NULL
  )
}

#' Prepare the intervention-agnostic Track B working data
#'
#' @description
#' Build the lag columns (`lag1_A`, `lag1_L`, ...) via
#' `causatr:::create_lag_vars()` and the risk-set bookkeeping
#' (`.survatr_prev_event` / `.survatr_prev_cens`) via `build_risk_set()` once
#' per `contrast()` call. Both are independent of the intervention and of the
#' horizon, so the per-(intervention, horizon) backward passes reuse them.
#'
#' @param fit A `survatr_fit` with `track == "B"`.
#' @param details Output of `build_ice_surv_details()`.
#'
#' @return A list with `data_lag` (person-period `data.table` with lag +
#'   `.survatr_prev_*` columns, keyed by `(id, time)`), `fit_rows` (logical
#'   at-risk mask aligned to `data_lag`), and `cols` (column-name bundle).
#' @noRd
prepare_track_b_base <- function(fit, details) {
  pp <- data.table::copy(fit$pp_data)
  data.table::setkeyv(pp, c(fit$id, fit$time))

  ## Lag every treatment + TV-confounder column. `create_lag_vars()` deep-
  ## copies and re-keys, so `data_lag` is a fresh object safe to mutate with
  ## the risk-set columns below.
  data_lag <- causatr:::create_lag_vars(
    pp,
    fit$treatment,
    details$tv_vars,
    fit$id,
    fit$time,
    fit$history
  )

  ## At-risk mask: prev_event == 0 (not yet failed), prev_cens == 0 and
  ## uncensored at this period. This mask defines BOTH the per-step fit set
  ## and -- crucially for the variance reuse -- the at-risk `fit_ids` that
  ## make causatr's forward sensitivity chain valid for the survival-tail
  ## pseudo-outcome (the `(1 - D_k)` factor is absorbed by restricting to
  ## the at-risk set; see `.claude/hard-rules.md`).
  fit_rows <- build_risk_set(
    data = data_lag,
    outcome = fit$outcome,
    id = fit$id,
    censoring = fit$censoring
  )

  list(
    data_lag = data_lag,
    fit_rows = fit_rows,
    cols = list(
      id = fit$id,
      time = fit$time,
      outcome = fit$outcome,
      treatment = fit$treatment,
      censoring = fit$censoring
    )
  )
}

#' Run one ICE backward pass to a fixed horizon
#'
#' @description
#' Estimate the counterfactual cumulative risk `R^d(τ)` for a single horizon
#' `τ = time_points[horizon_pos]` by Robins' iterated conditional expectation
#' on the **hazard link**, extended to time-to-event via the survival-tail
#' pseudo-outcome. Reuses causatr's per-step primitives wholesale; survatr
#' supplies only the survival target, the per-step link forcing, and the
#' at-risk fit sets.
#'
#' @details
#' Backward sweep from period `τ` (the fixed final step) down to the first
#' period:
#'
#' 1. **Final step (`p == horizon_pos`).** Fit the all-cause hazard among the
#'    at-risk-at-`τ` rows on the **observed** event indicator (binomial),
#'    then predict the counterfactual `q_{i,τ}` for those rows under the
#'    intervention (`data_iv`: current-period treatment intervened, lag
#'    columns left observed per `causatr:::ice_apply_intervention_long`).
#' 2. **Backward step (`p < horizon_pos`).** Form the survival-tail
#'    pseudo-outcome `Ỹ_k = ice_surv_pseudo(D_k, q_next)` among at-risk-at-`k`
#'    rows, fit it (quasibinomial), and predict the counterfactual `q_{i,k}`.
#'    Survived-but-unpredictable rows (`q_next` is `NA`, e.g. censored at the
#'    next period with no IPCW -- deferred to chunk 11) are dropped from the
#'    fit.
#' 3. After the sweep, `pseudo` at the first period equals the individual
#'    counterfactual `Ỹ_{0,i}(τ) = 1 - Ŝ^d_{0:τ,i}`; its mean is `R^d(τ)`.
#'
#' The returned list is shape-compatible with `causatr:::ice_iterate()`'s
#' return so it can be fed directly to `causatr:::variance_if_ice_one()`.
#'
#' @param data_lag Person-period data with lag + `.survatr_prev_*` columns
#'   (from `prepare_track_b_base()`).
#' @param data_iv Intervention-modified copy of `data_lag` (current-period
#'   treatment intervened, lags observed) from
#'   `causatr:::ice_apply_intervention_long()`. Reused across horizons of the
#'   same intervention.
#' @param fit_rows Logical at-risk mask aligned to `data_lag`.
#' @param horizon_pos Integer position of the horizon in `details$time_points`.
#' @param details Output of `build_ice_surv_details()`.
#' @param cols Column-name bundle (`id`, `time`, `outcome`, `treatment`,
#'   `censoring`).
#' @param intervention The `causatr_intervention` (carried into the result for
#'   the variance bridge).
#'
#' @return A 5-field list (`pseudo_final`, `models`, `fit_ids`, `data_iv`,
#'   `intervention`) -- the `ice_result` contract consumed by
#'   `causatr:::variance_if_ice_one()`.
#' @noRd
run_ice_survival_horizon <- function(
  data_lag,
  data_iv,
  fit_rows,
  horizon_pos,
  details,
  cols,
  intervention
) {
  id_col <- cols$id
  time_col <- cols$time
  outcome <- cols$outcome
  treatment <- cols$treatment

  time_points <- details$time_points
  n_times <- details$n_times
  model_fn <- details$model_fn
  dots <- details$dots
  external_weights <- details$weights

  time_vec <- data_lag[[time_col]]
  id_vec_all <- as.character(data_lag[[id_col]])

  ## Rolling per-individual pseudo-outcome, keyed by id. Filled at the final
  ## step and overwritten at each earlier step with the counterfactual tail
  ## prediction `q_{i,k}`.
  all_ids <- unique(id_vec_all)
  pseudo <- stats::setNames(rep(NA_real_, length(all_ids)), all_ids)

  ## Per-step model + fit-id storage, indexed by step (1 = earliest period)
  ## so the indexing matches `details$time_points` and the variance chain.
  ## Steps after the horizon stay NULL; the chain skips them.
  models <- vector("list", n_times)
  names(models) <- as.character(time_points)
  fit_ids <- vector("list", n_times)
  names(fit_ids) <- as.character(time_points)

  for (p in seq(horizon_pos, 1L, by = -1L)) {
    t_p <- time_points[p]
    time_idx <- p - 1L # 0-based for ice_build_formula's lag accounting

    ## At-risk-at-p rows: the fit set AND the prediction set.
    fit_mask <- fit_rows & (time_vec == t_p)
    fit_data <- data.table::copy(data_lag[fit_mask])
    d_p <- fit_data[[outcome]]
    ids_p <- as.character(fit_data[[id_col]])
    keep <- rep(TRUE, length(ids_p))

    if (p == horizon_pos) {
      ## Final step: response is the observed 0/1 event indicator; binomial.
      response_col <- outcome
      family_k <- details$family_outcome
    } else {
      ## Backward step: survival-tail pseudo-outcome; quasibinomial. `q_next`
      ## is the carried prediction from the later step (set in `pseudo`).
      q_next <- pseudo[ids_p]
      y_tilde <- ice_surv_pseudo(d_p, q_next)
      ## Drop survived-but-unpredictable rows (NA pseudo: at risk at p but not
      ## at p+1, i.e. censored next period with no IPCW). Failed-at-p rows
      ## keep their deterministic 1.
      keep <- !is.na(y_tilde)
      fit_data <- fit_data[keep]
      ids_p <- ids_p[keep]
      fit_data[, .pseudo_y := y_tilde[keep]]
      response_col <- ".pseudo_y"
      family_k <- details$family_pseudo
    }

    ## Per-step formula with lag + EM expansion (reused from causatr).
    formula_k <- causatr:::ice_build_formula(
      response_col,
      treatment,
      details$baseline_terms,
      details$tv_vars,
      details$tv_terms,
      time_idx,
      details$max_lag,
      fit_data,
      details$em_info,
      details$treatment_terms
    )

    ## External weights (Track B v1: always NULL) subset to the kept fit rows.
    w_k <- if (!is.null(external_weights)) {
      external_weights[which(fit_mask)][keep]
    } else {
      NULL
    }

    ## Fit the per-step model (pooled; `stratified = NULL`). `ice_fit_step()`
    ## centralises the family / weights / dots handling.
    models[[p]] <- causatr:::ice_fit_step(
      model_fn,
      formula_k,
      fit_data,
      family_k,
      w_k,
      dots,
      NULL
    )
    ## `fit_ids[[p]]` MUST be exactly the rows model `p` was fit on, in fit
    ## order -- the variance chain aligns per-step scores to ids through this.
    fit_ids[[p]] <- ids_p

    ## Predict the counterfactual tail `q_{i,p}` for ALL at-risk-at-p rows
    ## (the standardisation population at step p), conditioning on observed
    ## history with the current-period treatment intervened. Predict on the
    ## full at-risk set (pre-`keep`) so a row dropped from THIS fit still gets
    ## a pseudo for the earlier step.
    pred_data <- data_iv[fit_mask]
    preds <- causatr:::ice_predict_step(models[[p]], pred_data, NULL)
    pred_ids <- as.character(pred_data[[id_col]])
    pseudo[pred_ids] <- preds
  }

  ## `pseudo_final`: first-period ids in first-period row order (everyone is
  ## at risk at baseline on a rectangular grid), so `mean()` is `R^d(τ)`.
  first_t <- time_points[1]
  first_rows <- which(time_vec == first_t)
  first_ids <- id_vec_all[first_rows]
  pseudo_final <- unname(pseudo[first_ids])

  list(
    pseudo_final = pseudo_final,
    models = models,
    data_iv = data_iv,
    fit_ids = fit_ids,
    intervention = intervention
  )
}

#' Compute the Track B survival curve for one intervention
#'
#' @description
#' Driver over the requested horizons: one independent ICE backward pass per
#' `t` (the survival-tail pseudo-outcome is defined relative to a fixed final
#' horizon, so each `t` is its own terminal step). Returns the estimates table
#' (same shape as Track A's `compute_survival_curve()`) plus the per-horizon
#' `ice_result`s for the sandwich bridge.
#'
#' @param base Output of `prepare_track_b_base()`.
#' @param details Output of `build_ice_surv_details()`.
#' @param intervention A single `causatr_intervention`.
#' @param times User horizon grid (sorted unique, validated upstream; all in
#'   `details$time_points`).
#' @param intervention_name Name carried into the `intervention` column.
#'
#' @return A list with `estimates` (`data.table` with columns
#'   `intervention | time | s_hat | risk_hat | se | ci_lower | ci_upper | n`,
#'   `se` / `ci_*` as `NA_real_`) and `ice_results` (per-horizon list).
#' @noRd
compute_ice_survival_curve <- function(
  base,
  details,
  intervention,
  times,
  intervention_name
) {
  cols <- base$cols
  ## Build the intervention-modified frame ONCE (horizon-independent): current-
  ## period treatment intervened, lag columns left at observed values.
  data_iv <- causatr:::ice_apply_intervention_long(
    base$data_lag,
    cols$treatment,
    intervention,
    cols$id,
    cols$time
  )

  ## Standardisation population = individuals AT RISK at the first period.
  ## Individuals censored at entry (period 1) are never in the period-1 risk
  ## set, never receive a period-1 prediction, and so carry `NA` in
  ## `pseudo_final`. Under (MCAR) entry censoring they drop from the ICE
  ## standardisation -- the consistent g-formula behaviour -- so the mean is
  ## taken with `na.rm` over the at-risk-at-baseline ids. (Without a censoring
  ## column every first-period id is at risk, so `n_eff == n_ids` and this is a
  ## no-op.) Issue #1, 2026-06-03 critical review (/tmp/survatr_repro_cens.R).
  first_t <- details$time_points[1]
  first_mask <- base$data_lag[[cols$time]] == first_t
  n_eff <- sum(base$fit_rows[first_mask])

  ice_results <- vector("list", length(times))
  risk <- numeric(length(times))
  for (j in seq_along(times)) {
    horizon_pos <- match(times[j], details$time_points)
    res <- run_ice_survival_horizon(
      data_lag = base$data_lag,
      data_iv = data_iv,
      fit_rows = base$fit_rows,
      horizon_pos = horizon_pos,
      details = details,
      cols = cols,
      intervention = intervention
    )
    ice_results[[j]] <- res
    risk[j] <- mean(res$pseudo_final, na.rm = TRUE)
  }

  estimates <- data.table::data.table(
    intervention = intervention_name,
    time = times,
    s_hat = 1 - risk,
    risk_hat = risk,
    se = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    n = n_eff
  )

  list(estimates = estimates, ice_results = ice_results)
}
