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
#' @param type One of `"survival"`, `"risk"`, `"risk_difference"`,
#'   `"risk_ratio"`, `"rmst"`, `"rmst_difference"`.
#' @param reference Name of the intervention used as the reference in
#'   difference / ratio contrasts. Defaults to the first name in
#'   `interventions`. Ignored by `type = "survival"`, `"risk"`, and
#'   `"rmst"` (no pairwise contrast).
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
  type = c(
    "risk_difference",
    "survival",
    "risk",
    "risk_ratio",
    "rmst",
    "rmst_difference"
  ),
  reference = NULL,
  ci_method = "none",
  conf_level = 0.95,
  n_boot = 500L,
  boot_ci = "percentile",
  parallel = "no",
  ncpus = 1L,
  seed = NULL,
  ...
) {
  type <- match.arg(type)

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
  if (
    type %in%
      c("risk_difference", "risk_ratio", "rmst_difference") &&
      length(interventions) < 2L
  ) {
    rlang::abort(
      paste0(
        "`type = \"",
        type,
        "\"` requires at least two interventions (one reference + one ",
        "comparator). Pass a second intervention, or use a curve-only ",
        "type like \"survival\", \"risk\", or \"rmst\"."
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

  ## If the user asked for RMST-shaped output, replace the per-time s_hat
  ## column with the cumulative trapezoidal integral of S from 0 to t.
  if (type %in% c("rmst", "rmst_difference")) {
    estimates <- add_rmst_to_estimates(estimates, times)
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
        ipw_corr = shared$ipw_corr
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
      n_ids = length(shared$unique_ids)
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
    call = match.call()
  )
}
