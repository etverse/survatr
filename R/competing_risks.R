#' Fit parallel cause-specific hazard models (competing risks)
#'
#' Internal fit-only helper for the competing-risks path. For `J` competing
#' event types it fits `J` pooled-logistic cause-specific hazard models
#' `logit h^(j)(t | A, L) = alpha_j(t) + beta_{A,j} A + beta_{L,j} L` on a
#' **shared** all-cause risk set.
#'
#' The shared risk set is the key to "treat other causes as censored at their
#' event time": `build_risk_set()` is called upstream on the all-cause event
#' indicator `1{competing != 0}`, so an individual leaves the risk set at the
#' first event of **any** cause. Each cause-`j` model then regresses the
#' cause-`j` indicator `1{competing == j}` on those at-risk rows -- a row at a
#' competing cause `j' != j` event is still in the risk set (the event happens
#' *at* that period) with cause-`j` outcome 0, and is dropped only from the next
#' period on. This reproduces the cause-specific hazard exactly.
#'
#' @param data Full person-period `data.table` (after risk-set construction;
#'   includes the `.survatr_*` bookkeeping columns). Mutated in place with a
#'   transient `.survatr_cause_event` column, reset per cause; stripped later.
#' @param fit_rows Logical mask from `build_risk_set()` on the all-cause event
#'   indicator (shared across all `J` cause models).
#' @param competing Column name (character scalar) of the multi-valued cause
#'   column (`0` = no event, `1..J` = cause).
#' @param treatment,confounders Treatment column name and the baseline
#'   confounder RHS formula, used identically for every cause model.
#' @param time_formula RHS formula for the baseline hazard `alpha_j(t)`, shared
#'   across causes.
#' @param model_fn Fitting function (e.g. `stats::glm`, `mgcv::gam`).
#' @param causes Integer vector of the positive cause labels (from
#'   `check_competing_col()`).
#' @param ... Forwarded to `model_fn` (checked upstream to reject `na.exclude`).
#'
#' @returns A list with `model` (`NULL` -- a competing-risks fit has no single
#'   hazard model), `cause_models` (named list of the `J` fitted models, names
#'   `as.character(causes)`), `causes`, `family_name` (`"binomial"`), and
#'   `n_fit` (the shared at-risk row count).
#' @family competing-risks
#' @noRd
fit_competing_risks <- function(
  data,
  fit_rows,
  competing,
  treatment,
  confounders,
  time_formula,
  model_fn = stats::glm,
  causes,
  ...
) {
  cause_models <- vector("list", length(causes))
  names(cause_models) <- as.character(causes)
  for (j in causes) {
    ## Cause-`j` event indicator. Reset in place each iteration; each glm
    ## captures its own model frame at fit time, so reuse is safe.
    data[, .survatr_cause_event := as.integer(get(competing) == j)]
    fj <- fit_hazard_gcomp(
      data = data,
      fit_rows = fit_rows,
      outcome = ".survatr_cause_event",
      treatment = treatment,
      confounders = confounders,
      time_formula = time_formula,
      weights = NULL,
      model_fn = model_fn,
      ...
    )
    cause_models[[as.character(j)]] <- fj$model
  }

  list(
    model = NULL,
    cause_models = cause_models,
    causes = causes,
    family_name = "binomial",
    n_fit = sum(fit_rows)
  )
}

#' Attach per-cause hazards, all-cause hazard, and (lagged) survival to PP data
#'
#' Predicts each cause-specific hazard on the counterfactual person-period data,
#' attaches them as `.cf_haz_<j>` columns, then derives the all-cause hazard
#' `.cf_all_haz = sum_j h^(j)`, the within-id cumulative survival
#' `.cf_surv = prod_{m<=k} (1 - H_m)`, and the lagged survival
#' `.cf_surv_lag = S(k-1)` (with `S(0) = 1`). Centralised here so the survival
#' curve, every cause's CIF, and the influence-function chain all read the same
#' columns in one canonical `(id, time)` order.
#'
#' @param pp_cf Counterfactual person-period `data.table` (from
#'   `apply_intervention_pp()`). Mutated and returned.
#' @param cause_models Named list of fitted cause-specific hazard models.
#' @param causes Integer vector of cause labels.
#' @param id,time Column names for the id and time variables.
#'
#' @returns `pp_cf`, keyed by `(id, time)`, with the `.cf_haz_<j>`,
#'   `.cf_all_haz`, `.cf_surv`, and `.cf_surv_lag` columns attached.
#' @family competing-risks
#' @noRd
cr_augment_pp <- function(pp_cf, cause_models, causes, id, time) {
  data.table::setkeyv(pp_cf, c(id, time))
  haz_cols <- paste0(".cf_haz_", causes)
  for (j in causes) {
    hj <- predict_hazard_pp(cause_models[[as.character(j)]], pp_cf)
    pp_cf[, (paste0(".cf_haz_", j)) := hj]
  }
  ## All-cause hazard is the sum of cause-specific hazards (never average
  ## hazards before cumulating -- Jensen). `.SD` scoped to the per-cause hazard
  ## columns keeps the row-wise sum vectorised.
  pp_cf[, .cf_all_haz := rowSums(.SD), .SDcols = haz_cols]
  ## Within-id cumulative survival S_i(k) and its one-period lag S_i(k-1).
  pp_cf[, .cf_surv := cumprod(1 - .cf_all_haz), by = c(id)]
  pp_cf[,
    .cf_surv_lag := data.table::shift(.cf_surv, n = 1L, fill = 1, type = "lag"),
    by = c(id)
  ]
  pp_cf
}

#' Per-cause cumulative incidence curves under one intervention
#'
#' For each requested cause `j`, build the per-individual CIF
#' `F^(j)_i(t) = sum_{k<=t} S_i(k-1) h^(j)_{i,k}` (cumulative within id) and
#' average across ids at each requested time. Operates on the augmented PP data
#' from `cr_augment_pp()`.
#'
#' @param pp_cf Augmented counterfactual PP `data.table` (`.cf_haz_<j>`,
#'   `.cf_surv_lag` columns present).
#' @param id,time Column names.
#' @param times Numeric vector of requested time points (validated upstream).
#' @param causes Integer vector of the causes to report.
#' @param intervention_name Character scalar carried into the `intervention`
#'   column.
#'
#' @returns A `data.table` with columns `intervention | cause | time | cif_hat |
#'   se | ci_lower | ci_upper | n`, one row per (cause, time). `se` / `ci_*` are
#'   `NA_real_` placeholders filled by the variance path.
#' @family competing-risks
#' @noRd
compute_cif_curves <- function(
  pp_cf,
  id,
  time,
  times,
  causes,
  intervention_name
) {
  time_col <- time
  pieces <- lapply(causes, function(j) {
    hcol <- paste0(".cf_haz_", j)
    ## Per-id cumulative incidence increment dF = S(k-1) * h^(j), accumulated
    ## within id. `get(hcol)` pulls the cause-`j` hazard column inside the
    ## grouped expression (dynamic name; no bare-symbol global needed).
    pp_cf[, .cf_cif := cumsum(.cf_surv_lag * get(hcol)), by = c(id)]
    agg <- pp_cf[
      get(time_col) %in% times,
      list(cif_hat = mean(.cf_cif), n = .N),
      by = c(time_col)
    ]
    data.table::setnames(agg, time_col, "time")
    agg[,
      `:=`(
        intervention = intervention_name,
        cause = j,
        se = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
    ]
    agg
  })
  res <- data.table::rbindlist(pieces)
  data.table::setcolorder(
    res,
    c(
      "intervention",
      "cause",
      "time",
      "cif_hat",
      "se",
      "ci_lower",
      "ci_upper",
      "n"
    )
  )
  res[]
}

#' All-cause survival curve from the summed cause-specific hazards
#'
#' Averages the within-id all-cause survival `S_i(t) = prod_{m<=t} (1 - H_m)`
#' (with `H = sum_j h^(j)`) across ids at each requested time. Carries a
#' `cause = NA` marker so a competing-risks `survatr_result` keeps a stable
#' `cause` column across all its rows.
#'
#' @param pp_cf Augmented counterfactual PP `data.table` (`.cf_surv` present).
#' @param id,time Column names (`id` unused beyond keying upstream; kept for a
#'   symmetric signature with `compute_cif_curves()`).
#' @param times Numeric vector of requested time points.
#' @param intervention_name Character scalar.
#'
#' @returns A `data.table` with columns `intervention | cause | time | s_hat |
#'   risk_hat | se | ci_lower | ci_upper | n` (one row per time; `cause` is
#'   `NA_integer_`).
#' @family competing-risks
#' @noRd
compute_cr_survival_curve <- function(
  pp_cf,
  id,
  time,
  times,
  intervention_name
) {
  time_col <- time
  agg <- pp_cf[
    get(time_col) %in% times,
    list(s_hat = mean(.cf_surv), risk_hat = 1 - mean(.cf_surv), n = .N),
    by = c(time_col)
  ]
  data.table::setnames(agg, time_col, "time")
  agg[,
    `:=`(
      intervention = intervention_name,
      cause = NA_integer_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  ]
  data.table::setcolorder(
    agg,
    c(
      "intervention",
      "cause",
      "time",
      "s_hat",
      "risk_hat",
      "se",
      "ci_lower",
      "ci_upper",
      "n"
    )
  )
  agg[]
}

#' Assemble the competing-risks `contrasts` table
#'
#' `cif_difference` / `cif_ratio` compare each non-reference intervention's CIF
#' against the reference at every (cause, time). Curve-only types (`cif`,
#' `survival`, `risk`) return an empty stub with the canonical CR columns so
#' downstream S3 methods see a stable shape.
#'
#' @param estimates Per-intervention CIF `data.table` (with a `cause` column).
#' @param type One of the competing-risks contrast types.
#' @param reference Reference intervention name, or `NULL`.
#' @param interventions Named list of interventions.
#' @param causes Integer vector of the reported causes.
#'
#' @returns A `data.table` with columns `contrast | cause | time | estimate |
#'   se | ci_lower | ci_upper`. Empty (0 rows) for curve-only types.
#' @family competing-risks
#' @noRd
build_cif_contrasts <- function(
  estimates,
  type,
  reference,
  interventions,
  causes
) {
  empty_contrasts <- data.table::data.table(
    contrast = character(0),
    cause = integer(0),
    time = numeric(0),
    estimate = numeric(0),
    se = numeric(0),
    ci_lower = numeric(0),
    ci_upper = numeric(0)
  )
  if (estimand_is_curve(type)) {
    return(empty_contrasts)
  }

  other_names <- setdiff(names(interventions), reference)
  if (length(other_names) == 0L) {
    return(empty_contrasts)
  }

  ref_slim <- estimates[
    get("intervention") == reference,
    c("cause", "time", "cif_hat"),
    with = FALSE
  ]
  data.table::setnames(ref_slim, "cif_hat", "ref_val")

  pieces <- lapply(other_names, function(a1_name) {
    a1_slim <- estimates[
      get("intervention") == a1_name,
      c("cause", "time", "cif_hat"),
      with = FALSE
    ]
    data.table::setnames(a1_slim, "cif_hat", "a1_val")
    ## Join on (cause, time): each is unique within an intervention.
    merged <- merge(ref_slim, a1_slim, by = c("cause", "time"))
    est <- switch(
      type,
      cif_difference = merged$a1_val - merged$ref_val,
      cif_ratio = merged$a1_val / merged$ref_val
    )
    data.table::data.table(
      contrast = paste0(a1_name, " vs ", reference),
      cause = merged$cause,
      time = merged$time,
      estimate = est,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  })
  data.table::rbindlist(pieces)
}

#' Emit the competing-risks truncation-by-death caveat (once per session)
#'
#' Cause-specific CIF contrasts carry a well-known interpretational hazard: a
#' difference / ratio in cause-`j` cumulative incidence conditions on surviving
#' the competing events, so it is not a contrast of the same target population
#' under a sharp-null competing process. We surface this with a one-time session
#' message (rather than per call, to avoid spam) and the `print` method repeats
#' it for difference / ratio results. Returns the caveat text so callers can
#' attach it where useful.
#'
#' @returns The caveat string, invisibly emitted via `rlang::inform()` once per
#'   session.
#' @family competing-risks
#' @noRd
cr_truncation_caveat <- function() {
  msg <- paste0(
    "Cause-specific CIF contrasts condition on surviving the competing ",
    "events (truncation by death). Interpret a per-cause CIF difference / ",
    "ratio as a total effect on cause-specific incidence, not as a contrast ",
    "holding the competing process fixed. See the competing-risks vignette."
  )
  rlang::inform(
    msg,
    class = "survatr_cr_truncation_caveat",
    .frequency = "once",
    .frequency_id = "survatr_cr_truncation_caveat"
  )
  invisible(msg)
}

#' Competing-risks contrast path
#'
#' Orchestrates the cause-specific hazards + CIF contrast: builds the
#' counterfactual PP data per intervention, attaches per-cause hazards and the
#' all-cause survival via `cr_augment_pp()`, computes either the per-cause CIF
#' curves (`cif` / `cif_difference` / `cif_ratio`) or the all-cause survival /
#' risk curve (`survival` / `risk`), assembles contrasts, and fills the
#' sandwich or bootstrap variance. Mirrors `contrast.survatr_fit()` for the
#' single-event path but carries a `cause` dimension throughout.
#'
#' @param fit A competing-risks `survatr_fit` (`cause_models` non-`NULL`).
#' @param interventions Named list of `causatr_intervention` objects.
#' @param times Validated numeric time grid.
#' @param type One of `"cif"`, `"cif_difference"`, `"cif_ratio"`, `"survival"`,
#'   `"risk"`.
#' @param reference Reference intervention name, or `NULL` (curve-only types).
#' @param cause Integer vector of causes to report, or `NULL` (all). Validated
#'   against `fit$causes`; ignored for `survival` / `risk`.
#' @param ci_method,conf_level,n_boot,boot_ci,parallel,ncpus,seed Variance
#'   controls, forwarded from `contrast()`.
#' @param call The `match.call()` of `contrast()` for the result object.
#'
#' @returns A `survatr_result` with `cause`-bearing `estimates` and `contrasts`.
#' @family competing-risks
#' @noRd
contrast_competing <- function(
  fit,
  interventions,
  times,
  type,
  reference,
  cause,
  q_vec,
  ci_method,
  conf_level,
  n_boot,
  boot_ci,
  parallel,
  ncpus,
  seed,
  call
) {
  causes <- validate_cause(cause, fit$causes, type)
  ## CIF estimands carry the per-cause dimension; survival / risk are all-cause.
  is_cif <- estimand_has_cause(type)

  estimates_list <- lapply(names(interventions), function(iv_name) {
    iv <- interventions[[iv_name]]
    pp_cf <- apply_intervention_pp(fit$pp_data, fit$treatment, iv)
    pp_cf <- cr_augment_pp(
      pp_cf,
      fit$cause_models,
      fit$causes,
      fit$id,
      fit$time
    )
    if (is_cif) {
      compute_cif_curves(pp_cf, fit$id, fit$time, times, causes, iv_name)
    } else {
      compute_cr_survival_curve(pp_cf, fit$id, fit$time, times, iv_name)
    }
  })
  estimates <- data.table::rbindlist(estimates_list)

  ## All-cause survival quantile on a competing-risks fit. `type = "quantile"`
  ## carries no cause dimension (like survival / risk), so `is_cif` is FALSE and
  ## `estimates` already holds the all-cause `s_hat`. Build the all-cause CR
  ## survival IF (under sandwich) and hand off to the shared quantile assembly.
  if (identical(type, "quantile")) {
    iv_names <- names(interventions)
    s_by_iv <- survival_curves_by_iv(estimates, iv_names)
    if_list <- NULL
    n_ids <- length(unique(fit$pp_data[[fit$id]]))
    if (identical(ci_method, "sandwich")) {
      shared <- prepare_cr_sandwich_shared(fit)
      n_ids <- length(shared$unique_ids)
      if_list <- lapply(iv_names, function(iv) {
        pieces <- cr_intervention_if_pieces(fit, interventions[[iv]], shared)
        compute_cr_survival_if_matrix(pieces, times)$IF_mat
      })
      names(if_list) <- iv_names
    }
    return(finalize_quantile(
      fit = fit,
      interventions = interventions,
      times = times,
      q_vec = q_vec,
      reference = reference,
      s_by_iv = s_by_iv,
      if_list = if_list,
      n_ids = n_ids,
      ci_method = ci_method,
      conf_level = conf_level,
      n_boot = n_boot,
      boot_ci = boot_ci,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed,
      call = call
    ))
  }

  ## All-cause RMST / RMTL on a competing-risks fit: the same time-indexed
  ## transforms of the all-cause survival curve as the single-event path (`is_cif`
  ## is FALSE for these, so `estimates` holds `s_hat`).
  if (type %in% c("rmst")) {
    estimates <- add_rmst_to_estimates(estimates, times)
  }
  if (type %in% c("rmtl")) {
    estimates <- add_rmtl_to_estimates(estimates, times)
  }
  ## Years of life lost: integrate each cause's CIF (`is_cif` TRUE, so `estimates`
  ## holds per-cause `cif_hat`).
  if (identical(type, "yll")) {
    estimates <- add_yll_to_estimates(estimates, times)
  }

  contrasts <- build_cif_contrasts(
    estimates = estimates,
    type = type,
    reference = reference,
    interventions = interventions,
    causes = causes
  )

  ## Truncation-by-death caveat for the conditional CIF contrasts.
  if (estimand_has_cause(type) && estimand_is_contrast(type)) {
    cr_truncation_caveat()
  }

  if (identical(ci_method, "sandwich")) {
    shared <- prepare_cr_sandwich_shared(fit)
    filled <- fill_sandwich_ses_cr(
      fit = fit,
      estimates = estimates,
      contrasts = contrasts,
      interventions = interventions,
      times = times,
      type = type,
      reference = reference,
      causes = causes,
      conf_level = conf_level,
      shared = shared
    )
    estimates <- filled$estimates
    contrasts <- filled$contrasts
  } else if (identical(ci_method, "bootstrap")) {
    boot_out <- bootstrap_survival(
      fit = fit,
      interventions = interventions,
      times = times,
      type = type,
      reference = reference,
      n_boot = n_boot,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed,
      causes = causes
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
