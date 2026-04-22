#' Contrast generic
#'
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
#' @export
contrast <- function(fit, ...) {
  UseMethod("contrast")
}

#' Survival-curve contrasts on a fitted hazard model
#'
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
  ## Pairwise contrast types need at least two interventions (a non-
  ## reference vs a reference). Reject the single-intervention case
  ## upfront with a clear signal rather than silently returning an empty
  ## contrasts table or erroring deep in the replicate pipeline.
  if (
    type %in% c("risk_difference", "risk_ratio", "rmst_difference") &&
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
        unique_ids = shared$unique_ids
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

#' Validate the `interventions` argument
#'
#' Must be a non-empty named list; every element must inherit from
#' `causatr_intervention`; names must be unique and syntactically usable.
#'
#' @param interventions Value passed to `contrast()`.
#' @param call Caller frame for the error signal.
#'
#' @return Invisibly `NULL`.
#' @noRd
validate_interventions <- function(interventions, call = rlang::caller_env()) {
  if (!is.list(interventions) || length(interventions) == 0L) {
    rlang::abort(
      "`interventions` must be a non-empty named list.",
      class = "survatr_bad_interventions",
      call = call
    )
  }
  nms <- names(interventions)
  if (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    rlang::abort(
      "`interventions` must have unique, non-empty names.",
      class = "survatr_bad_interventions",
      call = call
    )
  }
  is_iv <- vapply(
    interventions,
    inherits,
    logical(1L),
    what = "causatr_intervention"
  )
  if (!all(is_iv)) {
    rlang::abort(
      c(
        "Every element of `interventions` must be a `causatr_intervention`.",
        i = paste0(
          "Construct with `causatr::static()`, `causatr::shift()`, ",
          "`causatr::scale_by()`, `causatr::threshold()`, or ",
          "`causatr::dynamic()`."
        )
      ),
      class = "survatr_bad_interventions",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the `times` argument
#'
#' Must be a numeric vector, no NAs, all values in `fit$time_grid`. Sorts
#' ascending and drops duplicates. Extrapolation beyond the observed time
#' grid is rejected -- chunk 2 does not handle it. Later chunks may relax
#' this to a warning with a dedicated classed signal.
#'
#' @param times User-supplied time grid.
#' @param fit_grid `fit$time_grid` (sorted unique observed times).
#' @param call Caller frame.
#'
#' @return The sorted, de-duplicated `times` vector.
#' @noRd
validate_times <- function(times, fit_grid, call = rlang::caller_env()) {
  ## Allow numeric / integer time grids (the common case) as well as
  ## Date / POSIXct / difftime, which users get when the hazard model
  ## uses a real-world timestamp column. `setdiff()` works across these
  ## types, so delegation to `fit_grid` for set membership is safe once
  ## the structural checks (non-empty, no NA) pass.
  is_time_like <- is.numeric(times) ||
    inherits(times, c("Date", "POSIXct", "POSIXlt", "difftime"))
  if (!is_time_like || length(times) == 0L || anyNA(times)) {
    rlang::abort(
      paste0(
        "`times` must be a non-empty vector of numeric / Date / POSIXct / ",
        "difftime values with no NA entries."
      ),
      class = "survatr_bad_times",
      call = call
    )
  }
  times <- sort(unique(times))
  missing <- setdiff(times, fit_grid)
  if (length(missing) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`times` contains value(s) not present in `fit$time_grid`: ",
          paste(missing, collapse = ", "),
          "."
        ),
        i = paste0(
          "survatr chunk 2 does not extrapolate beyond observed periods. ",
          "Pass only values in `fit$time_grid` = [",
          min(fit_grid),
          ", ",
          max(fit_grid),
          "]."
        )
      ),
      class = "survatr_time_extrapolation",
      call = call
    )
  }
  times
}

#' Validate the `reference` argument
#'
#' For difference / ratio contrasts, the reference must name one of the
#' interventions. For curve-only contrast types (`survival`, `risk`,
#' `rmst`), `reference` is ignored and returned as `NULL`.
#'
#' @param reference User-supplied reference name or `NULL`.
#' @param interventions The named list of interventions.
#' @param type Contrast type.
#' @param call Caller frame.
#'
#' @return Validated reference name, or `NULL` for curve-only types.
#' @noRd
validate_reference <- function(
  reference,
  interventions,
  type,
  call = rlang::caller_env()
) {
  no_contrast <- type %in% c("survival", "risk", "rmst")
  if (no_contrast) {
    return(NULL)
  }
  nms <- names(interventions)
  if (is.null(reference)) {
    reference <- nms[1L]
  }
  if (!reference %in% nms) {
    rlang::abort(
      paste0(
        "`reference = \"",
        reference,
        "\"` is not a name in `interventions` (",
        paste0("\"", nms, "\"", collapse = ", "),
        ")."
      ),
      class = "survatr_bad_reference",
      call = call
    )
  }
  reference
}

#' Validate the `ci_method` argument
#'
#' Accepts `"none"`, `"sandwich"`, and `"bootstrap"`. Anything else is
#' `survatr_bad_ci_method`.
#'
#' @param ci_method Scalar character.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @noRd
validate_ci_method <- function(ci_method, call = rlang::caller_env()) {
  if (ci_method %in% c("none", "sandwich", "bootstrap")) {
    return(invisible(NULL))
  }
  rlang::abort(
    paste0(
      "`ci_method` must be one of \"none\", \"sandwich\", \"bootstrap\". Got \"",
      ci_method,
      "\"."
    ),
    class = "survatr_bad_ci_method",
    call = call
  )
}

#' Validate the `n_boot` argument
#'
#' Must be a positive integer (coerced from numeric scalars when the
#' user passes e.g. `500`).
#'
#' @noRd
validate_n_boot <- function(n_boot, call = rlang::caller_env()) {
  if (
    !is.numeric(n_boot) ||
      length(n_boot) != 1L ||
      is.na(n_boot) ||
      n_boot < 1 ||
      abs(n_boot - round(n_boot)) > 1e-8
  ) {
    rlang::abort(
      paste0(
        "`n_boot` must be a positive integer. Got ",
        deparse(n_boot),
        "."
      ),
      class = "survatr_bad_n_boot",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the `boot_ci` argument
#'
#' @noRd
validate_boot_ci <- function(boot_ci, call = rlang::caller_env()) {
  if (!boot_ci %in% c("percentile", "wald")) {
    rlang::abort(
      paste0(
        "`boot_ci` must be one of \"percentile\", \"wald\". Got \"",
        boot_ci,
        "\"."
      ),
      class = "survatr_bad_boot_ci",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate `parallel` / `ncpus`
#'
#' Forwards to `boot::boot()`. Must be one of `"no"`, `"multicore"`,
#' `"snow"`; ncpus must be a positive integer.
#'
#' @noRd
validate_parallel <- function(parallel, ncpus, call = rlang::caller_env()) {
  if (!parallel %in% c("no", "multicore", "snow")) {
    rlang::abort(
      paste0(
        "`parallel` must be one of \"no\", \"multicore\", \"snow\". Got \"",
        parallel,
        "\"."
      ),
      class = "survatr_bad_parallel",
      call = call
    )
  }
  if (
    !is.numeric(ncpus) ||
      length(ncpus) != 1L ||
      is.na(ncpus) ||
      ncpus < 1 ||
      abs(ncpus - round(ncpus)) > 1e-8
  ) {
    rlang::abort(
      paste0("`ncpus` must be a positive integer. Got ", deparse(ncpus), "."),
      class = "survatr_bad_parallel",
      call = call
    )
  }
  invisible(NULL)
}

#' Validate the `conf_level` argument
#'
#' @param conf_level Numeric scalar.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @noRd
validate_conf_level <- function(conf_level, call = rlang::caller_env()) {
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1L ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    rlang::abort(
      paste0(
        "`conf_level` must be a single numeric value in (0, 1). Got ",
        deparse(conf_level),
        "."
      ),
      class = "survatr_bad_conf_level",
      call = call
    )
  }
  invisible(NULL)
}
