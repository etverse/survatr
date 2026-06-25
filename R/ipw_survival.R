#' Fit the IPW weighted hazard MSM (point-treatment, ipw)
#'
#' @description
#' Internal fit-only helper for `estimator = "ipw"`. Estimates the
#' counterfactual survival curve by inverse-probability weighting rather than
#' outcome standardization: fit a baseline treatment model, form stabilized
#' density-ratio weights at the observed treatment, broadcast the per-id weight
#' onto every person-period row, and fit a **weighted, marginal**
#' pooled-logistic hazard MSM `logit h(t | A) = alpha(t) + beta_A A` (no
#' confounders -- the confounding is handled by the weights).
#'
#' @details
#' The stabilized weight for individual `i` is
#'
#' ```
#' w_i = f(A_i) / f(A_i | L_i)
#' ```
#'
#' the ratio of the marginal treatment density to the conditional treatment
#' density, both evaluated at the **observed** treatment `A_i` (Hernán & Robins
#' 2020, Ch. 17). The weight is intervention-free and constant within `id`, so
#' it is broadcast unchanged onto every person-period row of `i`. The MSM is fit
#' with `quasibinomial()` (the weighted-fit family switch), via the existing
#' `fit_hazard_gcomp()` helper with `confounders = ~ 1` so the right-hand side
#' collapses to `alpha(t) + A`.
#'
#' causatr is the engine: the treatment model, the densities, and the
#' winsorization are all causatr internals (`fit_treatment_model()`,
#' `evaluate_density()`, `truncate_weights()`). survatr only composes the
#' stabilized ratio and broadcasts it onto the person-period grid. causatr's own
#' point IPW uses Hájek normalization on a `Y ~ 1` MSM and never forms this
#' observed-treatment stabilized weight, so there is no causatr function to
#' reuse for the composition itself.
#'
#' This chunk supports **binary** treatment only; continuous / categorical /
#' count propensities are a dedicated future path and abort here with
#' `survatr_ipw_treatment_unsupported`.
#'
#' @param data Full person-period `data.table` (validated, sorted by
#'   `(id, time)`, includes the `.survatr_*` bookkeeping columns).
#' @param fit_rows Logical mask from `build_risk_set()` selecting the at-risk
#'   rows for the hazard MSM.
#' @param outcome,treatment Column names (character scalars).
#' @param confounders One-sided formula for the **treatment** model
#'   (`A ~ confounders`). The hazard MSM uses `A` only.
#' @param id,time Column names of the id and discrete-period index.
#' @param time_formula One-sided formula for the baseline hazard `alpha(t)`.
#' @param propensity_model_fn Fitting function for the treatment model
#'   (`stats::glm`-compatible). Default `stats::glm`.
#' @param trim Numeric scalar in `(0, 1]` or `NULL`. When `< 1`, the per-id
#'   weights are winsorized at the `trim`-th quantile (Cole & Hernán 2008). The
#'   resolved fixed cutoff is returned so the sandwich uses the same threshold.
#' @param model_fn Fitting function for the hazard MSM. Default `stats::glm`.
#' @param ... Forwarded to `model_fn` for the hazard MSM (same convention as the
#'   gcomp path). Not forwarded to the treatment model.
#'
#' @returns A list with `model` (the fitted weighted hazard MSM),
#'   `family_name`, `n_fit`, `treatment_model` (the full `A ~ L` causatr
#'   treatment model), `marginal_model` (the `A ~ 1` numerator model), `weights`
#'   (the per-PP-row broadcast weight vector, length `nrow(data)`), and
#'   `trim_threshold` (the fixed winsorization cutoff, or `NA_real_`).
#' @family survatr_fit functions
#' @noRd
fit_ipw_survival <- function(
  data,
  fit_rows,
  outcome,
  treatment,
  confounders,
  id,
  time,
  time_formula,
  propensity_model_fn = stats::glm,
  trim = NULL,
  model_fn = stats::glm,
  ...
) {
  ## Collapse to one baseline row per id. `data` is keyed by (id, time), so the
  ## first row within each id is the earliest period -- the baseline at which
  ## the point treatment and its confounders are measured. `.I[1L]` returns the
  ## global row index of that first row per id; subsetting keeps the id sort.
  baseline_idx <- data[, .I[1L], by = c(id)]$V1
  baseline <- data[baseline_idx]

  ## The point treatment must be constant within id for point-treatment
  ## g-computation: the weight is a single per-id quantity broadcast across
  ## periods. A time-varying treatment here means the user has longitudinal data
  ## (a longitudinal-IPW problem), which this path cannot weight correctly.
  trt_uniq <- data[,
    list(.nuniq = data.table::uniqueN(.SD[[1L]])),
    by = c(id),
    .SDcols = treatment
  ]
  if (any(trt_uniq$.nuniq > 1L)) {
    rlang::abort(
      c(
        paste0(
          "IPW (`estimator = \"ipw\"`) requires a point treatment that is ",
          "constant within `",
          id,
          "`, but `",
          treatment,
          "` varies within at least one id."
        ),
        i = paste0(
          "Time-varying treatment is a longitudinal-IPW problem and is not ",
          "handled by the point-treatment hazard MSM."
        )
      ),
      class = "survatr_ipw_time_varying_treatment"
    )
  }

  ## Positivity / no-variation guard. If every id shares the same treatment
  ## value there is no contrast to weight toward, and causatr's family
  ## auto-detection would see a single unique value and misclassify the binary
  ## treatment as "gaussian" (its bernoulli check needs exactly two values),
  ## producing a misleading "continuous treatment" error. Catch it here with a
  ## clear positivity message instead.
  if (data.table::uniqueN(baseline[[treatment]]) < 2L) {
    rlang::abort(
      c(
        paste0(
          "IPW requires variation in `",
          treatment,
          "`, but every individual shares the same treatment value."
        ),
        i = paste0(
          "With no treated/untreated contrast the inverse-probability ",
          "weights are undefined (a positivity violation)."
        )
      ),
      class = "survatr_ipw_no_treatment_variation"
    )
  }

  ## Full propensity model A ~ L on the baseline rows. causatr fits the density,
  ## drops aliased (collinear) columns, and stores the design `X_prop` and the
  ## aligned `alpha_hat` we reuse in the weight closure for the sandwich.
  full_tm <- causatr:::fit_treatment_model(
    data = baseline,
    treatment = treatment,
    confounders = confounders,
    model_fn = propensity_model_fn
  )

  ## Binary-only gate for this chunk. Continuous (gaussian) / categorical /
  ## count propensities need pushforward / multinomial weight forms and their
  ## own oracles; they ship in a dedicated extended-types path.
  if (!identical(full_tm$family, "bernoulli")) {
    rlang::abort(
      c(
        paste0(
          "IPW survival currently supports a binary treatment only; `",
          treatment,
          "` was detected as family \"",
          full_tm$family,
          "\"."
        ),
        i = paste0(
          "Continuous / categorical / count treatment IPW for survival ",
          "ships in a dedicated extended-types path."
        )
      ),
      class = "survatr_ipw_treatment_unsupported"
    )
  }

  ## Marginal numerator model A ~ 1 (stabilization). For a binary treatment this
  ## is the marginal P(A = 1); fitting it via the same causatr machinery keeps
  ## the density evaluation uniform and aligned with the denominator model.
  marginal_model <- causatr:::fit_treatment_model(
    data = baseline,
    treatment = treatment,
    confounders = ~1,
    model_fn = propensity_model_fn
  )

  ## Stabilized per-id weights at the observed treatment, winsorized at a fixed
  ## cutoff. The cutoff is returned so the sandwich winsorizes at the same
  ## threshold (the bootstrap re-estimates it per replicate via a fresh fit).
  wres <- compute_ipw_stabilized_weights(
    full_tm = full_tm,
    marginal_model = marginal_model,
    baseline = baseline,
    treatment = treatment,
    trim = trim
  )

  ## Broadcast the per-id weight onto every person-period row by id.
  weights_pp <- broadcast_weights_to_pp(
    w_by_id = wres$weights,
    baseline_ids = baseline[[id]],
    pp_ids = data[[id]]
  )

  ## Weighted marginal hazard MSM: alpha(t) + A, quasibinomial (weighted fit).
  ## `confounders = ~ 1` collapses the RHS to the time basis plus the treatment.
  haz <- fit_hazard_gcomp(
    data = data,
    fit_rows = fit_rows,
    outcome = outcome,
    treatment = treatment,
    confounders = ~1,
    time_formula = time_formula,
    weights = weights_pp,
    model_fn = model_fn,
    ...
  )

  list(
    model = haz$model,
    family_name = haz$family_name,
    n_fit = haz$n_fit,
    treatment_model = full_tm,
    marginal_model = marginal_model,
    weights = weights_pp,
    trim_threshold = wres$trim_threshold,
    ## IPCW fields: NULL for the standard IPW-only path.
    censoring_model = NULL,
    ipcw_numerator_model = NULL,
    ipcw_weights = NULL,
    ipcw_trim_thresholds = NULL,
    ipw_treatment_weights_pp = NULL
  )
}

#' Fit the IPW + IPCW weighted hazard MSM (point-treatment, ipw + built-in IPCW)
#'
#' @description
#' Extends `fit_ipw_survival()` with a per-period inverse-probability-of-
#' censoring weight. When censoring depends on measured covariates, the
#' pure row-filter path produces a biased curve; the IPCW weight
#' `W^C_{i,k} = prod_{m <= k} P(C_m = 0 | A) / P(C_m = 0 | A, L, m)`
#' up-weights survivors at each period so the weighted hazard MSM estimates
#' the uncensored potential-outcome curve.
#'
#' The procedure:
#' 1. Estimate the stabilized treatment weight `w_i = f(A_i) / f(A_i | L_i)`
#'    (same as `fit_ipw_survival()`).
#' 2. Fit the censoring hazard model `P(C_k = 1 | A, L, time)` on the
#'    at-risk, not-yet-censored rows.
#' 3. Compute the per-period cumulative IPCW weight `W^C_{i,k}` as a running
#'    product.
#' 4. Combine: row weight for PP row `(i, k)` = `w_i * W^C_{i,k}`.
#' 5. Fit the weighted marginal hazard MSM with the combined weight.
#'
#' The two weight components are stored separately in the returned list so the
#' sandwich variance can fix the treatment component while perturbing the
#' censoring-model parameters (and vice versa).
#'
#' @param data Full person-period `data.table` (validated, sorted, with
#'   `.survatr_prev_*` columns already written by `build_risk_set()`).
#' @param fit_rows Logical mask from `build_risk_set()`.
#' @param outcome,treatment,censoring Column names.
#' @param confounders One-sided formula for the treatment model.
#' @param id,time Column names.
#' @param time_formula One-sided time-basis formula.
#' @param propensity_model_fn Fitting function for the treatment model.
#' @param trim Quantile for per-id treatment weight winsorization AND
#'   per-period IPCW weight winsorization.
#' @param model_fn Fitting function for the hazard MSM.
#' @param ipcw_formula One-sided formula for the censoring-model covariates.
#' @param censoring_model_fn Fitting function for the censoring hazard.
#' @param ... Forwarded to `model_fn` and `censoring_model_fn`.
#'
#' @returns A list with `model`, `family_name`, `n_fit`, `treatment_model`,
#'   `marginal_model`, `weights` (combined per-PP-row weight),
#'   `trim_threshold`, `censoring_model`, `ipcw_numerator_model`,
#'   `ipcw_weights` (per-PP-row cumulative censoring weights),
#'   `ipcw_trim_thresholds`, and `ipw_treatment_weights_pp`
#'   (treatment-only broadcast weights).
#' @family survatr_fit functions
#' @noRd
fit_ipw_survival_ipcw <- function(
  data,
  fit_rows,
  outcome,
  treatment,
  confounders,
  id,
  time,
  time_formula,
  propensity_model_fn,
  trim,
  model_fn,
  ipcw_formula,
  censoring_model_fn,
  censoring,
  ...
) {
  ## Step 1: treatment weights (same as fit_ipw_survival()).
  baseline_idx <- data[, .I[1L], by = c(id)]$V1
  baseline <- data[baseline_idx]

  trt_uniq <- data[,
    list(.nuniq = data.table::uniqueN(.SD[[1L]])),
    by = c(id),
    .SDcols = treatment
  ]
  if (any(trt_uniq$.nuniq > 1L)) {
    rlang::abort(
      c(
        paste0(
          "IPW + IPCW requires a point treatment that is constant within `",
          id,
          "`, but `",
          treatment,
          "` varies within at least one id."
        ),
        i = "Time-varying treatment is a longitudinal-IPW problem."
      ),
      class = "survatr_ipw_time_varying_treatment"
    )
  }
  if (data.table::uniqueN(baseline[[treatment]]) < 2L) {
    rlang::abort(
      c(
        paste0(
          "IPW requires variation in `",
          treatment,
          "`, but every individual shares the same treatment value."
        ),
        i = paste0(
          "With no treated/untreated contrast the inverse-probability ",
          "weights are undefined (a positivity violation)."
        )
      ),
      class = "survatr_ipw_no_treatment_variation"
    )
  }

  full_tm <- causatr:::fit_treatment_model(
    data = baseline,
    treatment = treatment,
    confounders = confounders,
    model_fn = propensity_model_fn
  )
  if (!identical(full_tm$family, "bernoulli")) {
    rlang::abort(
      c(
        paste0(
          "IPW + IPCW currently supports a binary treatment only; `",
          treatment,
          "` was detected as family \"",
          full_tm$family,
          "\"."
        ),
        i = paste0(
          "Continuous / categorical / count treatment IPW ships in a ",
          "dedicated extended-types path."
        )
      ),
      class = "survatr_ipw_treatment_unsupported"
    )
  }
  marginal_model <- causatr:::fit_treatment_model(
    data = baseline,
    treatment = treatment,
    confounders = ~1,
    model_fn = propensity_model_fn
  )
  wres <- compute_ipw_stabilized_weights(
    full_tm = full_tm,
    marginal_model = marginal_model,
    baseline = baseline,
    treatment = treatment,
    trim = trim
  )
  ## Broadcast per-id treatment weight onto all PP rows (constant within id).
  ipw_weights_pp <- broadcast_weights_to_pp(
    w_by_id = wres$weights,
    baseline_ids = baseline[[id]],
    pp_ids = data[[id]]
  )

  ## Step 2: fit the censoring hazard model and compute per-period IPCW weights.
  ## `build_risk_set()` has already been called by surv_fit() and written
  ## `.survatr_prev_*` columns onto `data` in-place, so we can call
  ## `fit_censoring_model()` directly.
  cres <- fit_censoring_model(
    data = data,
    censoring = censoring,
    treatment = treatment,
    ipcw_formula = ipcw_formula,
    time_formula = time_formula,
    censoring_model_fn = censoring_model_fn,
    id = id,
    time = time,
    ...
  )
  wc_res <- compute_ipcw_running_weights(
    data = data,
    cens_model = cres$cens_model,
    num_model = cres$num_model,
    id = id,
    time = time,
    trim = trim
  )

  ## Step 3: combine treatment weight × IPCW weight for each PP row.
  combined_weights_pp <- ipw_weights_pp * wc_res$weights

  ## Step 4: fit the weighted marginal hazard MSM with the combined weight.
  haz <- fit_hazard_gcomp(
    data = data,
    fit_rows = fit_rows,
    outcome = outcome,
    treatment = treatment,
    confounders = ~1,
    time_formula = time_formula,
    weights = combined_weights_pp,
    model_fn = model_fn,
    ...
  )

  list(
    model = haz$model,
    family_name = haz$family_name,
    n_fit = haz$n_fit,
    treatment_model = full_tm,
    marginal_model = marginal_model,
    ## Combined weight stored in `weights` (used by the hazard MSM) so
    ## existing contrast / variance machinery reads it unchanged.
    weights = combined_weights_pp,
    trim_threshold = wres$trim_threshold,
    ## IPCW-specific outputs stored separately for the censoring sandwich block.
    censoring_model = cres$cens_model,
    ipcw_numerator_model = cres$num_model,
    ipcw_weights = wc_res$weights,
    ipcw_trim_thresholds = wc_res$trim_thresholds,
    ipw_treatment_weights_pp = ipw_weights_pp
  )
}
