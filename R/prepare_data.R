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
  if (any(rows_per_id$.n == 1L)) {
    n_wide <- sum(rows_per_id$.n == 1L)
    rlang::abort(
      paste0(
        "Data does not appear to be in person-period format: ",
        n_wide,
        " of ",
        nrow(rows_per_id),
        " unique `",
        id,
        "` values have only a single row. ",
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

  data
}
