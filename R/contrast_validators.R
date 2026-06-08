#' Validate the `interventions` argument
#'
#' Must be a non-empty named list; every element must inherit from
#' `causatr_intervention`; names must be unique and syntactically usable.
#'
#' @param interventions Value passed to `contrast()`.
#' @param call Caller frame for the error signal.
#'
#' @return Invisibly `NULL`.
#' @family checks
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
#' @family checks
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
#' @family checks
#' @noRd
validate_reference <- function(
  reference,
  interventions,
  type,
  call = rlang::caller_env()
) {
  no_contrast <- type %in% c("survival", "risk", "rmst", "cif")
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

#' Validate the `cause` argument (competing risks)
#'
#' For CIF estimands (`cif` / `cif_difference` / `cif_ratio`), `cause` selects
#' which competing cause(s) to report. `NULL` means all fitted causes. Any
#' requested cause must be one of the cause labels present in the fit
#' (`fit$causes`). For the all-cause estimands (`survival` / `risk`) the `cause`
#' argument is meaningless and is ignored (returns `NULL`).
#'
#' @param cause User-supplied cause selection: `NULL` or an integer-valued
#'   vector.
#' @param fit_causes The fitted cause labels (`fit$causes`), or `NULL` for a
#'   single-event fit.
#' @param type Contrast type.
#' @param call Caller frame for the error signal.
#'
#' @returns Integer vector of the causes to report (sorted, de-duplicated), or
#'   `NULL` for the all-cause estimands.
#' @family checks
#' @noRd
validate_cause <- function(
  cause,
  fit_causes,
  type,
  call = rlang::caller_env()
) {
  ## All-cause estimands carry no cause dimension.
  if (type %in% c("survival", "risk")) {
    return(NULL)
  }
  ## Default: report every fitted cause.
  if (is.null(cause)) {
    return(fit_causes)
  }
  if (!is.numeric(cause) || length(cause) == 0L || anyNA(cause)) {
    rlang::abort(
      "`cause` must be `NULL` or a non-empty integer vector with no NA.",
      class = "survatr_bad_cause",
      call = call
    )
  }
  cause <- sort(unique(as.integer(round(cause))))
  unknown <- setdiff(cause, fit_causes)
  if (length(unknown) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`cause` contains label(s) not among the fitted causes: ",
          paste(unknown, collapse = ", "),
          "."
        ),
        i = paste0(
          "Fitted causes are: ",
          paste(fit_causes, collapse = ", "),
          "."
        )
      ),
      class = "survatr_bad_cause",
      call = call
    )
  }
  cause
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
#' @family checks
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
#' @param n_boot Number of bootstrap replicates passed to `contrast()`.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @family checks
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
#' Accepts `"percentile"` and `"wald"`. Anything else is
#' `survatr_bad_boot_ci`.
#'
#' @param boot_ci Scalar character bootstrap CI flavour.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @family checks
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
#' @param parallel Scalar character backend passed to `boot::boot()`.
#' @param ncpus Number of CPUs for the parallel bootstrap.
#' @param call Caller frame.
#'
#' @return Invisibly `NULL`.
#' @family checks
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
#' @family checks
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
