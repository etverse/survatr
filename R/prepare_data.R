#' Validate and coerce person-period data
#'
#' Internal preflight for `surv_fit()`:
#'
#' 1. Coerce `data.frame` / `tibble` input to `data.table` via a copy, so we
#'    never mutate the caller's object in place.
#' 2. Check that every column referenced by the user (outcome, treatment, id,
#'    time, optional censoring) is actually present.
#' 3. Check that the data is in person-period (long) shape -- every `id` has at
#'    least two rows, and `(id, time)` pairs are unique. Wide or duplicated
#'    data produces a classed error pointing to `causatr::to_person_period()`.
#' 4. Sort by `(id, time)` so downstream within-id lagged cumsums are well
#'    defined.
#'
#' @param data User-supplied data.
#' @param outcome,treatment,id,time Column names.
#' @param censoring Column name, or `NULL`.
#' @param call Enclosing frame for the error signal.
#'
#' @return A `data.table` sorted by `(id, time)`. Not a view; safe to mutate.
#' @noRd
prepare_pp_data <- function(
  data,
  outcome,
  treatment,
  id,
  time,
  censoring = NULL,
  call = rlang::caller_env()
) {
  if (data.table::is.data.table(data)) {
    data <- data.table::copy(data)
  } else {
    data <- data.table::as.data.table(data)
  }

  required_cols <- c(outcome, treatment, id, time, censoring)
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0L) {
    rlang::abort(
      paste0(
        "Column(s) ",
        paste0("`", missing, "`", collapse = ", "),
        " not found in `data`."
      ),
      class = "survatr_col_not_found",
      call = call
    )
  }

  data.table::setkeyv(data, c(id, time))

  rows_per_id <- data[, list(.n = .N), by = c(id)]
  n_ids <- nrow(rows_per_id)
  n_times <- length(unique(data[[time]]))

  ## Reject wide (one row per id) input. A degenerate edge case where
  ## the overall time grid has length 1 AND every id has 1 row is
  ## well-formed (trivially rectangular) and falls through; the more
  ## common "one row per id with k > 1 time values" is the wide shape
  ## and gets the pointer to `causatr::to_person_period()`.
  if (n_times > 1L && any(rows_per_id$.n == 1L)) {
    n_wide <- sum(rows_per_id$.n == 1L)
    rlang::abort(
      paste0(
        "Data does not appear to be in person-period format: ",
        n_wide,
        " of ",
        n_ids,
        " unique `",
        id,
        "` values have only a single row, and the overall time grid has ",
        n_times,
        " periods. ",
        "Use `causatr::to_person_period()` to convert wide data to long ",
        "before calling `surv_fit()`."
      ),
      class = "survatr_not_person_period",
      call = call
    )
  }

  ## Each id must have distinct time points. Duplicated (id, time) pairs
  ## would make the within-id lagged cumsums ill-defined (the lag would skip
  ## duplicates instead of respecting the period index).
  time_unique_per_id <- data[,
    list(nt = length(unique(.SD[[1L]]))),
    by = c(id),
    .SDcols = time
  ]
  if (any(time_unique_per_id$nt < rows_per_id$.n)) {
    rlang::abort(
      paste0(
        "Some individuals have duplicated rows at the same `",
        time,
        "` value. `surv_fit()` requires unique (",
        id,
        ", ",
        time,
        ") pairs."
      ),
      class = "survatr_duplicate_pp_row",
      call = call
    )
  }

  ## surv_fit() requires a rectangular PP grid: every id must have a row
  ## at every unique time value. The contrast / sandwich / bootstrap
  ## paths all pull one row per id at each requested time
  ## (`compute_survival_curve`, `compute_survival_if_matrix`,
  ## `bootstrap_survival`). Ragged PP (ids dropped post-event or
  ## post-censor) must be padded upfront: typically by appending rows
  ## with the outcome set to 0 and the censoring column (if any) set to
  ## 1, so the risk-set builder drops them from the fit while the
  ## prediction path still gets hazards at every (id, t) pair.
  if (nrow(data) != n_ids * n_times) {
    rlang::abort(
      c(
        paste0(
          "Ragged person-period data: ",
          nrow(data),
          " rows across ",
          n_ids,
          " ids and ",
          n_times,
          " unique times (expected ",
          n_ids * n_times,
          ")."
        ),
        i = paste0(
          "`surv_fit()` requires a rectangular PP grid. Pad each id to ",
          "the full time grid before calling `surv_fit()`; set the ",
          "outcome to 0 and the censoring column (if any) to 1 on the ",
          "padded rows so the risk-set builder drops them from the fit."
        )
      ),
      class = "survatr_ragged_pp",
      call = call
    )
  }

  data
}
