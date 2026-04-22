#' Build the pooled-logistic risk set
#'
#' Adds within-id lagged cumulative sums of the event (and censoring) indicators
#' directly onto `data` and returns the logical vector selecting the rows that
#' enter the hazard GLM fit.
#'
#' The risk set at period `k` for individual `i` is defined by:
#'
#' - `.survatr_prev_event_{i,k} = 0` -- the individual has not yet experienced
#'   the event by the end of period `k-1`.
#' - (if `censoring` non-NULL) `.survatr_prev_cens_{i,k} = 0` AND the censoring
#'   indicator at period `k` is `NA` or `0` (via `is_uncensored()`) -- the
#'   individual has not yet been censored and is not censored at this period.
#'
#' Lags are built with `data.table::shift(type = "lag", fill = 0)` after
#' sorting by `(id, time)`. Sorting must have happened upstream in
#' `prepare_pp_data()`.
#'
#' @param data Person-period `data.table`, already validated and sorted by
#'   `(id, time)`. Mutated in place with the `.survatr_prev_*` columns.
#' @param outcome,id,censoring Column names. `censoring` may be `NULL`.
#'
#' @return Logical vector of length `nrow(data)`. `TRUE` rows are the risk
#'   set -- the rows fed to the hazard GLM.
#' @noRd
build_risk_set <- function(data, outcome, id, censoring = NULL) {
  ## Lagged cumsum of events: at row (i, k), value is the number of events
  ## observed in periods < k. Zero means the individual is still at risk.
  data[,
    .survatr_prev_event := data.table::shift(
      cumsum(get(outcome)),
      n = 1L,
      fill = 0,
      type = "lag"
    ),
    by = c(id)
  ]

  if (!is.null(censoring)) {
    data[,
      .survatr_prev_cens := data.table::shift(
        cumsum(get(censoring)),
        n = 1L,
        fill = 0,
        type = "lag"
      ),
      by = c(id)
    ]
    fit_rows <- data[[".survatr_prev_event"]] == 0 &
      data[[".survatr_prev_cens"]] == 0 &
      is_uncensored(data, censoring)
  } else {
    fit_rows <- data[[".survatr_prev_event"]] == 0
  }

  fit_rows
}
