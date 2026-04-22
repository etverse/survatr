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
#' @param ci_method Only `"none"` is accepted at this chunk. `"sandwich"`
#'   and `"bootstrap"` ship in chunks 3 and 4; passing them now raises
#'   `survatr_ci_not_available`.
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
  ...
) {
  type <- match.arg(type)

  validate_interventions(interventions)
  times <- validate_times(times, fit$time_grid)
  reference <- validate_reference(reference, interventions, type)
  validate_ci_method(ci_method)

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
  if (!is.numeric(times) || length(times) == 0L || anyNA(times)) {
    rlang::abort(
      "`times` must be a non-empty numeric vector with no NA values.",
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

#' Reject unsupported ci_method values
#'
#' Only `"none"` is accepted at chunk 2. `"sandwich"` and `"bootstrap"` are
#' deliberately rejected so the user gets a signal pointing to the chunks
#' that wire them up rather than silently falling back to no CIs.
#'
#' @param ci_method Scalar character.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @noRd
validate_ci_method <- function(ci_method, call = rlang::caller_env()) {
  if (identical(ci_method, "none")) {
    return(invisible(NULL))
  }
  if (ci_method %in% c("sandwich", "bootstrap")) {
    rlang::abort(
      paste0(
        "`ci_method = \"",
        ci_method,
        "\"` is not available. survatr currently returns point estimates ",
        "only (`ci_method = \"none\"`)."
      ),
      class = "survatr_ci_not_available",
      call = call
    )
  }
  rlang::abort(
    paste0("`ci_method` must be \"none\". Got \"", ci_method, "\"."),
    class = "survatr_bad_ci_method",
    call = call
  )
}
