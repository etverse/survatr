#' Fit the IPW weighted hazard MSM (Track A, ipw)
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

  ## The point treatment must be constant within id for Track A: the weight is a
  ## single per-id quantity broadcast across periods. A time-varying treatment
  ## here means the user has longitudinal data (a Track B / longitudinal-IPW
  ## problem), which this path cannot weight correctly.
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
    trim_threshold = wres$trim_threshold
  )
}
