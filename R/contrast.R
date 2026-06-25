#' Contrast generic
#'
#' @description
#' S3 generic dispatched on the class of `fit`. survatr ships
#' `contrast.survatr_fit()` for time-indexed survival / risk / RMST
#' contrasts.
#'
#' `causatr::contrast()` (same symbol name in causatr, but not an S3
#' generic there) handles scalar-outcome causatr_fit objects. If both
#' packages are attached, this generic shadows `causatr::contrast()`;
#' users who need the causatr path should call `causatr::contrast()`
#' explicitly.
#'
#' @param fit A fitted model object.
#' @param ... Arguments passed to methods.
#'
#' @return Method-dependent.
#' @family survatr_result methods
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
#' fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
#'
#' # Contrast two static interventions; dispatches to contrast.survatr_fit().
#' contrast(
#'   fit,
#'   interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
#'   times = 1:4,
#'   type = "risk_difference"
#' )
#' @export
contrast <- function(fit, ...) {
  UseMethod("contrast")
}

#' Survival-curve contrasts on a fitted hazard model
#'
#' @description
#' Given a `survatr_fit` (chunk 1) and a named list of interventions, build
#' the counterfactual person-period data under each intervention, predict
#' the per-row hazard, cumulate within individual to get `S^a_i(t)`, and
#' average across individuals to get the population survival curve
#' `S^a(t) = (1/n) sum_i S^a_i(t)`. Derive risk / risk-difference /
#' risk-ratio / RMST / RMST-difference contrasts from there.
#'
#' **No variance in this chunk.** All `se` / `ci_*` columns come back as
#' `NA_real_`; sandwich variance ships with chunk 3 (delta-method cross-time
#' IF aggregation) and bootstrap with chunk 4.
#'
#' @param fit A `survatr_fit`.
#' @param interventions A named list of `causatr_intervention` objects
#'   (`causatr::static()`, `causatr::shift()`, `causatr::scale_by()`,
#'   `causatr::threshold()`, `causatr::dynamic()`). The element names are
#'   carried through to the `intervention` column of the result.
#' @param times Numeric vector of time points at which to evaluate the
#'   survival curve. Must all be elements of `fit$time_grid` — extrapolation
#'   beyond the observed grid is rejected with
#'   `survatr_time_extrapolation`.
#' @param type Estimand. For a single-event fit: one of `"survival"`,
#'   `"risk"`, `"risk_difference"`, `"risk_ratio"`, `"rmst"`,
#'   `"rmst_difference"`, `"rmtl"` (restricted mean time lost, `t - RMST`),
#'   `"rmtl_difference"`, or `"quantile"` (survival-time quantile / median;
#'   default `"risk_difference"`). For a competing-risks fit
#'   (`surv_fit(..., competing = )`): one of `"cif"` (per-cause cumulative
#'   incidence), `"cif_difference"`, `"cif_ratio"` (default), `"yll"` (per-cause
#'   years of life lost, `int F^(j)`), all-cause `"survival"` / `"risk"` /
#'   `"rmst"` / `"rmtl"`, or all-cause `"quantile"` (from the summed
#'   cause-specific hazards). Mixing a CIF estimand with a single-event fit (or
#'   vice versa) aborts with `survatr_competing_type`; `"yll"` on a single-event
#'   fit aborts with `survatr_yll_needs_cr`.
#' @param q For `type = "quantile"`, the quantile level(s) of the survival-time
#'   distribution to report -- a numeric vector with every entry in `(0, 1)`
#'   (default `0.5` = median). The result is indexed by `q` (one `tau_hat` per
#'   `(intervention, q)`) rather than by `time`. Aborts with
#'   `survatr_quantile_unreached` when the curve never crosses `1 - q` on the
#'   grid, and `survatr_bad_q` for out-of-range `q`. Ignored for other types.
#' @param cause Competing-risks only. Integer vector selecting which cause(s)
#'   to report for the `cif` estimands, or `NULL` (the default) for all causes.
#'   Validated against the fitted causes; ignored for `survival` / `risk` and
#'   for single-event fits.
#' @param reference Name of the intervention used as the reference in
#'   difference / ratio contrasts. Defaults to the first name in
#'   `interventions`. Ignored by `type = "survival"`, `"risk"`, `"rmst"`, and
#'   `"cif"` (no pairwise contrast).
#' @param ci_method One of `"none"` (point estimates only), `"sandwich"`
#'   (delta-method cross-time IF aggregation via
#'   `causatr:::prepare_model_if()`), or `"bootstrap"` (resample
#'   individuals, refit the hazard model per replicate, percentile or
#'   Wald CIs). Default `"none"`.
#' @param conf_level Confidence level for the CIs when `ci_method != "none"`.
#'   Numeric scalar in `(0, 1)`, default `0.95`.
#' @param n_boot Integer; number of bootstrap replicates. Ignored when
#'   `ci_method != "bootstrap"`. Default `500L`.
#' @param boot_ci One of `"percentile"` (sample-quantile CI) or `"wald"`
#'   (point estimate +/- `z * sd(replicates)`). Default `"percentile"`.
#'   Percentile is transform-invariant and is the safer default for
#'   ratios / RMST.
#' @param parallel One of `"no"`, `"multicore"`, `"snow"`; forwarded to
#'   `boot::boot()`. Default `"no"`.
#' @param ncpus Integer; number of CPUs for parallel bootstrap. Default
#'   `1L`.
#' @param seed Integer scalar or `NULL`. When non-null, `set.seed(seed)`
#'   is called before the bootstrap loop so the replicate sequence is
#'   reproducible. Default `NULL`.
#' @param cluster `NULL` (the default; per-individual sandwich) or a single
#'   column name in the fit's person-period data identifying the independent
#'   sampling cluster (site, household, provider, repeated enrolment). The label
#'   must be constant within `id` and have no NA. When supplied, the sandwich
#'   variance is cluster-robust (the per-individual influence functions are
#'   summed within cluster before the cross-product) and the bootstrap resamples
#'   whole clusters. `cluster = "<id-column>"` reproduces the per-individual
#'   sandwich exactly (singleton clusters). Supported for gcomp, IPW, IPCW,
#'   Track B (ICE), and competing-risks fits.
#' @param ... Unused; reserved for future chunks.
#'
#' @return A `survatr_result` list with `estimates`, `contrasts`,
#'   `time_grid`, `type`, `reference`, `ci_method`, and `call` components.
#' @family survatr_result methods
#' @examples
#' # Small rectangular person-period dataset: 40 ids over 5 periods.
#' set.seed(2)
#' n_id <- 40L
#' K <- 5L
#' pp <- data.frame(
#'   id = rep(seq_len(n_id), each = K),
#'   t = rep(seq_len(K), times = n_id),
#'   A = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   Y = rbinom(n_id * K, 1L, 0.1)
#' )
#'
#' fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
#'
#' ivs <- list(a1 = causatr::static(1), a0 = causatr::static(0))
#'
#' # Point-estimate survival curves under each intervention.
#' contrast(fit, interventions = ivs, times = 1:5, type = "survival")
#'
#' # Risk difference with delta-method sandwich CIs.
#' contrast(
#'   fit,
#'   interventions = ivs,
#'   times = 1:5,
#'   type = "risk_difference",
#'   ci_method = "sandwich"
#' )
#' @export
contrast.survatr_fit <- function(
  fit,
  interventions,
  times,
  type = NULL,
  cause = NULL,
  q = 0.5,
  reference = NULL,
  ci_method = "none",
  conf_level = 0.95,
  n_boot = 500L,
  boot_ci = "percentile",
  parallel = "no",
  ncpus = 1L,
  seed = NULL,
  cluster = NULL,
  ...
) {
  ## A competing-risks fit (J cause-specific hazards) and a single-event fit
  ## support disjoint estimand sets. Resolve the `type` default per fit kind --
  ## `cif_difference` for competing risks, `risk_difference` for single event --
  ## so `contrast(fit, ...)` with no `type` does the natural thing for each.
  is_competing <- !is.null(fit$cause_models)
  ## Valid-estimand allow-lists come from the central estimand registry
  ## (`estimand_registry.R`), so a new estimand is one descriptor entry rather
  ## than an edit to a hardcoded vector in every dispatch site.
  competing_types <- competing_estimands()
  if (is.null(type)) {
    type <- if (is_competing) "cif_difference" else "risk_difference"
  }
  type <- match.arg(type, estimand_types())

  ## The quantile estimand collapses the time axis onto a `q` dimension; `q` is
  ## only meaningful here. Validate elementwise so a vector `q = c(.25, .5,
  ## .75)` is accepted (one row per (intervention, q) in the result).
  q_vec <- if (identical(type, "quantile")) validate_q(q) else NULL

  ## Cross-check estimand vs fit kind. CIF estimands need a competing-risks
  ## fit; the single-event contrast estimands (risk_difference / risk_ratio /
  ## rmst*) are not defined per cause and are rejected for competing-risks fits
  ## (use cif_difference / cif_ratio, or all-cause survival / risk instead).
  if (identical(type, "yll") && !is_competing) {
    ## Years of life lost is a per-cause integral of the CIF, so it is only
    ## defined on a competing-risks fit. Distinct class from the generic CIF
    ## cross-check so callers can match the YLL-specific failure.
    rlang::abort(
      c(
        "`type = \"yll\"` (years of life lost) requires a competing-risks fit.",
        i = paste0(
          "Fit with `surv_fit(..., competing = <cause-col>)` first; ",
          "for a single-event fit use `type = \"rmtl\"`."
        )
      ),
      class = "survatr_yll_needs_cr"
    )
  }
  if (estimand_has_cause(type) && !is_competing) {
    rlang::abort(
      c(
        paste0("`type = \"", type, "\"` requires a competing-risks fit."),
        i = "Fit with `surv_fit(..., competing = <cause-col>)` first."
      ),
      class = "survatr_competing_type"
    )
  }
  if (is_competing && !type %in% competing_types) {
    rlang::abort(
      c(
        paste0(
          "`type = \"",
          type,
          "\"` is not defined for a competing-risks fit."
        ),
        i = paste0(
          "Use \"cif\", \"cif_difference\", or \"cif_ratio\" (per cause), or ",
          "\"survival\" / \"risk\" (all-cause)."
        )
      ),
      class = "survatr_competing_type"
    )
  }

  validate_interventions(interventions)
  ## IPSI (`causatr::ipsi()`) is an IPW-only intervention that reweights the
  ## propensity (Kennedy 2019) rather than plugging a fixed treatment value into
  ## the MSM, so it does not fit the static / dynamic plug-in path the survival
  ## curve is built on. Its weight-path survival estimand ships in a dedicated
  ## chunk; reject it here with a clear, classed message rather than letting the
  ## generic `apply_intervention()` abort surface.
  if (identical(fit$estimator, "ipw")) {
    iv_types <- vapply(
      interventions,
      function(iv) if (is.null(iv$type)) NA_character_ else iv$type,
      character(1L)
    )
    if (any(iv_types == "ipsi", na.rm = TRUE)) {
      rlang::abort(
        c(
          "`ipsi()` is not yet supported for IPW survival.",
          i = paste0(
            "Incremental propensity-score interventions use a weight-path ",
            "estimand (not an MSM plug-in) and ship in a dedicated chunk."
          )
        ),
        class = "survatr_ipw_ipsi_deferred"
      )
    }
  }
  ## Pairwise contrast types need at least two interventions (a non-
  ## reference vs a reference). Reject the single-intervention case
  ## upfront with a clear signal rather than silently returning an empty
  ## contrasts table or erroring deep in the replicate pipeline.
  if (estimand_requires_pair(type) && length(interventions) < 2L) {
    rlang::abort(
      paste0(
        "`type = \"",
        type,
        "\"` requires at least two interventions (one reference + one ",
        "comparator). Pass a second intervention, or use a curve-only ",
        "type like \"survival\", \"risk\", \"rmst\", or \"cif\"."
      ),
      class = "survatr_bad_interventions"
    )
  }
  times <- validate_times(times, fit$time_grid)
  reference <- validate_reference(reference, interventions, type)
  validate_ci_method(ci_method)
  validate_conf_level(conf_level)
  validate_n_boot(n_boot)
  validate_boot_ci(boot_ci)
  validate_parallel(parallel, ncpus)
  ## Cluster-robust variance: resolve the cluster column to per-id labels once
  ## (constant-within-id, no NA, G >= 2 are validated here). `NULL` keeps the
  ## per-individual sandwich. The labels are name-keyed (id -> cluster) so each
  ## variance path reindexes them onto its own IF row order.
  cluster_labels <- validate_cluster(cluster, fit)

  ## Competing risks (cause-specific hazards + CIF). Builds per-cause cumulative
  ## incidence (or all-cause survival) from the J cause-specific models and a
  ## stacked-EE sandwich across them. Carries a `cause` dimension the
  ## single-event curve path does not, so it takes a dedicated branch.
  if (is_competing) {
    return(contrast_competing(
      fit = fit,
      interventions = interventions,
      times = times,
      type = type,
      reference = reference,
      cause = cause,
      q_vec = q_vec,
      ci_method = ci_method,
      conf_level = conf_level,
      n_boot = n_boot,
      boot_ci = boot_ci,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed,
      cluster_labels = cluster_labels,
      call = match.call()
    ))
  }

  ## Track B (longitudinal ICE-hazard survival). The curve is built by a
  ## per-(intervention, horizon) backward sequential regression on the
  ## survival-tail pseudo-outcome rather than the Track A predict-hazard /
  ## cumulative-product path, so it takes a dedicated branch and returns early.
  ## The estimand shape, contrast assembly, RMST integral, and CI fillers are
  ## all reused -- only the curve + influence-function construction differ.
  if (identical(fit$track, "B")) {
    return(contrast_track_b(
      fit = fit,
      interventions = interventions,
      times = times,
      type = type,
      reference = reference,
      q_vec = q_vec,
      ci_method = ci_method,
      conf_level = conf_level,
      n_boot = n_boot,
      boot_ci = boot_ci,
      parallel = parallel,
      ncpus = ncpus,
      seed = seed,
      cluster_labels = cluster_labels,
      call = match.call()
    ))
  }

  ## Per-intervention survival curves. Build the counterfactual PP data,
  ## predict hazards on every row (at-risk rows are irrelevant for the
  ## prediction -- the cumulative product over k <= t needs the predicted
  ## hazard even at periods where the individual was censored / had the
  ## event in reality), compute the cumulative product within id, and
  ## average across ids at each t in `times`.
  estimates_list <- lapply(names(interventions), function(iv_name) {
    iv <- interventions[[iv_name]]
    pp_cf <- apply_intervention_pp(fit$pp_data, fit$treatment, iv)
    hazards <- predict_hazard_pp(fit$model, pp_cf)
    compute_survival_curve(
      pp_data = pp_cf,
      hazards = hazards,
      id = fit$id,
      time = fit$time,
      times = times,
      intervention_name = iv_name
    )
  })
  estimates <- data.table::rbindlist(estimates_list)

  ## Quantile (median / arbitrary q): a functional of the survival curve that
  ## collapses the time axis onto a `q` dimension, so it takes a dedicated
  ## branch and returns a q-indexed result. The IF matrices (under sandwich) are
  ## the same ones the time-indexed estimands use, so IPW / IPCW inherit the
  ## quantile through `compute_survival_if_matrix()` unchanged.
  if (identical(type, "quantile")) {
    iv_names <- names(interventions)
    s_by_iv <- survival_curves_by_iv(estimates, iv_names)
    if_list <- NULL
    cl_aligned <- NULL
    n_ids <- length(unique(fit$pp_data[[fit$id]]))
    if (identical(ci_method, "sandwich")) {
      shared <- prepare_sandwich_shared(fit)
      n_ids <- length(shared$unique_ids)
      if_list <- lapply(iv_names, function(iv_name) {
        compute_survival_if_matrix(
          fit = fit,
          intervention = interventions[[iv_name]],
          times = times,
          prep = shared$prep,
          fit_idx = shared$fit_idx,
          id_vec = shared$id_vec,
          unique_ids = shared$unique_ids,
          ipw_corr = shared$ipw_corr,
          ipcw_corr = shared$ipcw_corr
        )$IF_mat
      })
      names(if_list) <- iv_names
      ## Align the cluster labels onto the IF matrix row order for the sandwich.
      cl_aligned <- cluster_for_ids(cluster_labels, shared$unique_ids)
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
      call = match.call(),
      cluster_aligned = cl_aligned,
      cluster_labels = cluster_labels
    ))
  }

  ## If the user asked for RMST-shaped output, replace the per-time s_hat
  ## column with the cumulative trapezoidal integral of S from 0 to t.
  if (type %in% c("rmst", "rmst_difference")) {
    estimates <- add_rmst_to_estimates(estimates, times)
  }
  ## RMTL is the complement of RMST in the same window (t - RMST). Same
  ## time-indexed curve shape; the SE reuses the RMST quadratic form because
  ## the constant `t` drops out of the gradient.
  if (type %in% c("rmtl", "rmtl_difference")) {
    estimates <- add_rmtl_to_estimates(estimates, times)
  }

  ## Assemble pairwise contrasts (difference / ratio) at each t in `times`.
  ## For `survival` / `risk` / `rmst` the contrasts table is an empty stub
  ## with the canonical columns so S3 methods can dispatch on a stable shape.
  contrasts <- build_contrasts(
    estimates = estimates,
    type = type,
    reference = reference,
    interventions = interventions
  )

  ## Sandwich variance: delta-method cross-time IF aggregation. For each
  ## intervention we build the n_ids x |t| IF matrix on S^a(t), then
  ## propagate to the estimand's SE / CI via fill_sandwich_ses().
  if (identical(ci_method, "sandwich")) {
    shared <- prepare_sandwich_shared(fit)
    if_list <- lapply(names(interventions), function(iv_name) {
      compute_survival_if_matrix(
        fit = fit,
        intervention = interventions[[iv_name]],
        times = times,
        prep = shared$prep,
        fit_idx = shared$fit_idx,
        id_vec = shared$id_vec,
        unique_ids = shared$unique_ids,
        ipw_corr = shared$ipw_corr,
        ipcw_corr = shared$ipcw_corr
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
      n_ids = length(shared$unique_ids),
      ## Align the name-keyed cluster labels onto the IF matrix row order.
      cluster = cluster_for_ids(cluster_labels, shared$unique_ids)
    )
    estimates <- filled$estimates
    contrasts <- filled$contrasts
  } else if (identical(ci_method, "bootstrap")) {
    ## Empirical bootstrap: resample individuals, refit the hazard model
    ## per replicate, recompute curves / contrasts, percentile or Wald
    ## bands across replicates. Per-id resampling preserves the
    ## within-id cumulative-product dependence structure.
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
      cluster_labels = cluster_labels
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
    call = match.call()
  )
}
