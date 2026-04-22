## Column names survatr reserves for internal bookkeeping. Guarded upfront so
## user-data collisions produce an informative error rather than silent overwrite.
SURVATR_RESERVED_COLS <- c(".survatr_prev_event", ".survatr_prev_cens")

#' Validate external weights
#'
#' Copy-adapted from `causatr:::check_weights`. Accepts `NULL` (no weights).
#' Otherwise: must be numeric, length equal to `n`, finite, non-negative.
#' `NA` / `Inf` / `NaN` / negative values are rejected. Zero weights are
#' allowed (they act as row exclusions in GLM fits).
#'
#' @param weights Numeric vector, or `NULL`.
#' @param n Expected length (`nrow(data)` at the call site).
#' @param call Enclosing frame for the error signal; passed to `rlang::abort()`.
#'
#' @return Invisibly `NULL`. Called for its side effect (error or pass-through).
#' @noRd
check_weights <- function(weights, n, call = rlang::caller_env()) {
  if (is.null(weights)) {
    return(invisible(NULL))
  }
  if (!is.numeric(weights)) {
    rlang::abort(
      "`weights` must be numeric.",
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (length(weights) != n) {
    rlang::abort(
      paste0(
        "`weights` must have length equal to `nrow(data)` (",
        n,
        "), got ",
        length(weights),
        "."
      ),
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (anyNA(weights)) {
    rlang::abort(
      paste0(
        "`weights` contains ",
        sum(is.na(weights)),
        " missing value(s). ",
        "Drop those rows or impute before calling `surv_fit()`."
      ),
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (any(!is.finite(weights))) {
    rlang::abort(
      "`weights` contains non-finite value(s) (Inf / NaN).",
      class = "survatr_bad_weights",
      call = call
    )
  }
  if (any(weights < 0)) {
    rlang::abort(
      "`weights` must be non-negative.",
      class = "survatr_bad_weights",
      call = call
    )
  }
  invisible(NULL)
}

#' Reject `na.action = na.exclude` in `...`
#'
#' Copy-adapted from `causatr:::check_dots_na_action`. Under `na.exclude`,
#' `residuals(model, "working")` is padded with `NA`s while `model.matrix()`
#' silently drops them; the IF engine's `d_fit * r_score` product then
#' misaligns and silently corrupts the sandwich variance. We only accept
#' `na.omit` (default) or `na.fail`.
#'
#' @param ... Dots captured from the caller (e.g. `surv_fit()`).
#' @param call Enclosing frame for the error signal.
#'
#' @return Invisibly `NULL`. Called for its side effect.
#' @noRd
check_dots_na_action <- function(..., call = rlang::caller_env()) {
  dots <- list(...)
  if (!"na.action" %in% names(dots)) {
    return(invisible(NULL))
  }
  na_action <- dots$na.action
  ok <- FALSE
  if (is.function(na_action)) {
    ok <- identical(na_action, stats::na.omit) ||
      identical(na_action, stats::na.fail)
  } else if (is.character(na_action) && length(na_action) == 1L) {
    ok <- na_action %in% c("na.omit", "na.fail")
  }
  if (!ok) {
    rlang::abort(
      c(
        "`na.action` must be `na.omit` (default) or `na.fail`.",
        i = paste0(
          "survatr builds its own row-alignment bookkeeping from ",
          "`fit_rows` and the fitted model's `na.action` attribute. ",
          "`na.exclude` pads working residuals with NA and silently ",
          "corrupts the sandwich variance."
        ),
        i = "Drop NA rows before calling `surv_fit()` or use `na.action = na.omit`."
      ),
      class = "survatr_bad_na_action",
      call = call
    )
  }
  invisible(NULL)
}

#' Reject user-data collisions with survatr's reserved column names
#'
#' survatr adds `.survatr_prev_event` and `.survatr_prev_cens` to the
#' person-period data during risk-set construction. If the user's data
#' already contains a column with one of these names we would silently
#' overwrite it, so abort upfront.
#'
#' @param data Person-period data (anything with `names()`).
#' @param which Character vector of reserved names to check against.
#' @param call Enclosing frame for the error signal.
#'
#' @return Invisibly `NULL`. Called for its side effect.
#' @noRd
check_reserved_cols <- function(
  data,
  which = SURVATR_RESERVED_COLS,
  call = rlang::caller_env()
) {
  bad <- intersect(which, names(data))
  if (length(bad) > 0L) {
    rlang::abort(
      paste0(
        "Column name(s) ",
        paste0("`", bad, "`", collapse = ", "),
        " are reserved by survatr internals. Rename the column(s) in ",
        "your input data before calling `surv_fit()`."
      ),
      class = "survatr_reserved_col",
      call = call
    )
  }
  invisible(NULL)
}
