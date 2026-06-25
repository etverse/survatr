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
  ## Curve-only estimands carry no pairwise contrast, so `reference` is unused.
  if (estimand_is_curve(type)) {
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
  ## All-cause estimands (survival / risk) carry no cause dimension.
  if (!estimand_has_cause(type)) {
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

#' Validate the `q` argument (survival quantile level)
#'
#' For `type = "quantile"`, `q` selects which quantile(s) of the survival-time
#' distribution to report. Must be a non-empty numeric vector with every entry
#' strictly inside `(0, 1)` and no NA. A vector is accepted (one row per
#' `(intervention, q)` in the result).
#'
#' @param q User-supplied quantile level(s).
#' @param call Caller frame for the error signal.
#'
#' @returns The sorted, de-duplicated `q` vector.
#' @family checks
#' @noRd
validate_q <- function(q, call = rlang::caller_env()) {
  if (
    !is.numeric(q) ||
      length(q) == 0L ||
      anyNA(q) ||
      any(q <= 0) ||
      any(q >= 1)
  ) {
    rlang::abort(
      paste0(
        "`q` must be a non-empty numeric vector with every value strictly in ",
        "(0, 1). Got ",
        deparse(q),
        "."
      ),
      class = "survatr_bad_q",
      call = call
    )
  }
  sort(unique(q))
}

#' Validate the `cluster` argument and resolve per-id cluster labels
#'
#' @description
#' For cluster-robust variance, `cluster` names a column of the fit's
#' person-period data that labels the independent sampling unit (site,
#' household, provider, repeated enrolment). The cluster label must be constant
#' within `id` (an individual belongs to exactly one cluster) and have no NA;
#' there must be at least two clusters (a single cluster is not estimable). The
#' validated labels are collapsed to one value per `id` and returned as a
#' **named** character vector keyed by the id (as character), so each variance
#' path can reindex it onto its own influence-function row order via
#' `cluster_for_ids()`.
#'
#' @details
#' Returning the labels keyed by id (rather than aligned to a particular IF
#' matrix) decouples validation from row order: the single-event sandwich, the
#' competing-risks sandwich, the quantile assembly, and the bootstrap each have
#' their own id ordering, and they all reindex this one name-keyed vector. The
#' column is read but never mutated, so the fit's person-period data is left
#' untouched.
#'
#' @param cluster `NULL` (per-individual sandwich, the default) or a single
#'   column name in `fit$pp_data`.
#' @param fit A `survatr_fit` (supplies `pp_data` and the `id` column name).
#' @param call Caller frame for the error signal.
#'
#' @returns `NULL` when `cluster` is `NULL`, otherwise a named character vector
#'   of length `n_ids` mapping each id (as character) to its cluster label.
#' @family checks
#' @noRd
validate_cluster <- function(cluster, fit, call = rlang::caller_env()) {
  if (is.null(cluster)) {
    return(NULL)
  }
  ## A single existing column name. Anything else is a usage error rather than
  ## one of the data-shape failures below.
  if (
    !is.character(cluster) ||
      length(cluster) != 1L ||
      is.na(cluster) ||
      !nzchar(cluster)
  ) {
    rlang::abort(
      paste0(
        "`cluster` must be `NULL` or a single column name (character scalar). ",
        "Got ",
        deparse(cluster),
        "."
      ),
      class = "survatr_bad_cluster",
      call = call
    )
  }
  pp <- fit$pp_data
  if (!cluster %in% names(pp)) {
    rlang::abort(
      paste0(
        "`cluster = \"",
        cluster,
        "\"` is not a column in the fit's person-period data."
      ),
      class = "survatr_bad_cluster",
      call = call
    )
  }
  id_col <- fit$id
  ## NA cluster labels would silently drop ids from the within-cluster row sum
  ## (and `rowsum()` would create a phantom NA cluster), corrupting G and the
  ## meat -- reject upfront.
  if (anyNA(pp[[cluster]])) {
    rlang::abort(
      paste0(
        "`cluster = \"",
        cluster,
        "\"` has NA values; every row must carry a cluster label."
      ),
      class = "survatr_cluster_na",
      call = call
    )
  }
  ## An individual belongs to exactly one cluster. If a cluster label varies
  ## within an id, the per-id collapse below would be ambiguous and the
  ## row-sum unit-of-independence would be ill-defined.
  per_id_nu <- pp[,
    list(.nu = data.table::uniqueN(get(cluster))),
    by = c(id_col)
  ]
  if (any(per_id_nu$.nu > 1L)) {
    n_bad <- sum(per_id_nu$.nu > 1L)
    rlang::abort(
      c(
        paste0(
          "`cluster = \"",
          cluster,
          "\"` is not constant within `",
          id_col,
          "`: ",
          n_bad,
          " id(s) span more than one cluster."
        ),
        i = "Each individual must belong to exactly one cluster."
      ),
      class = "survatr_cluster_varies_within_id",
      call = call
    )
  }
  ## Collapse to one label per id, first-appearance order. Coerce to character
  ## so numeric / factor cluster columns key consistently and `rowsum()`'s
  ## factor() sees stable levels.
  per_id <- pp[,
    list(.cl = as.character(get(cluster))[1L]),
    by = c(id_col)
  ]
  labels <- per_id$.cl
  names(labels) <- as.character(per_id[[id_col]])
  ## Cluster-robust variance treats the cluster as the independent unit; with a
  ## single cluster there is exactly one unit and the between-cluster variance
  ## is undefined. Demand G >= 2.
  if (length(unique(labels)) < 2L) {
    rlang::abort(
      paste0(
        "`cluster = \"",
        cluster,
        "\"` resolves to a single cluster; cluster-robust variance needs at ",
        "least two clusters."
      ),
      class = "survatr_cluster_degenerate",
      call = call
    )
  }
  labels
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
      "`ci_method` must be one of \"none\", \"sandwich\", \"bootstrap\". ",
      "Got \"",
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
