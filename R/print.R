#' Print a `survatr_fit`
#'
#' Minimal banner summary for the fit object returned by `surv_fit()`.
#' Reports the track, estimator, outcome / treatment / id / time columns,
#' number of individuals, number of person-period rows used to fit, and the
#' time grid span. A richer print (coef table, time-spline degrees of
#' freedom) ships with the contrast / S3 polish in a later chunk.
#'
#' @param x A `survatr_fit`.
#' @param ... Unused.
#'
#' @return The fit object, invisibly.
#' @export
print.survatr_fit <- function(x, ...) {
  n_id <- length(unique(x$pp_data[[x$id]]))
  tg <- x$time_grid
  cens_info <- if (is.null(x$censoring)) "none" else x$censoring

  cat("<survatr_fit>\n")
  cat(sprintf("  Track:       %s\n", x$track))
  cat(sprintf("  Estimator:   %s\n", x$estimator))
  cat(sprintf("  Family:      %s\n", x$family))
  cat(sprintf("  Outcome:     %s\n", x$outcome))
  cat(sprintf("  Treatment:   %s\n", x$treatment))
  cat(sprintf("  ID:          %s\n", x$id))
  cat(sprintf("  Time:        %s\n", x$time))
  cat(sprintf("  Censoring:   %s\n", cens_info))
  cat(sprintf(
    "  N:           %d individuals, %d PP rows (%d at risk)\n",
    n_id,
    x$n_total,
    x$n_fit
  ))
  cat(sprintf(
    "  Time grid:   [%s, %s] (%d unique times)\n",
    format(tg[1L]),
    format(tg[length(tg)]),
    length(tg)
  ))
  invisible(x)
}

#' Print a `survatr_result`
#'
#' Minimal banner + head of the result's `contrasts` (or `estimates` for
#' curve-only `type`s). A polished print + `plot` / `tidy` / `forrest`
#' surface ships with the S3 polish in a later chunk.
#'
#' @param x A `survatr_result`.
#' @param n Maximum number of rows from the contrasts / estimates table to
#'   show (default 10).
#' @param ... Unused.
#'
#' @return The result object, invisibly.
#' @export
print.survatr_result <- function(x, n = 10L, ...) {
  tg <- x$time_grid
  cat("<survatr_result>\n")
  cat(sprintf("  Type:        %s\n", x$type))
  cat(sprintf(
    "  Reference:   %s\n",
    if (is.null(x$reference)) "(none)" else x$reference
  ))
  cat(sprintf("  CI method:   %s\n", x$ci_method))
  cat(sprintf(
    "  Time grid:   [%s, %s] (%d unique times)\n",
    format(tg[1L]),
    format(tg[length(tg)]),
    length(tg)
  ))
  cat(sprintf("  Estimates:   %d rows\n", nrow(x$estimates)))
  cat(sprintf("  Contrasts:   %d rows\n", nrow(x$contrasts)))

  show <- if (nrow(x$contrasts) > 0L) x$contrasts else x$estimates
  ## Use `head()` (not `show[seq_len(n)]`) because `show` may carry an
  ## `n` column (the per-time sample count from `compute_survival_curve`)
  ## that would bind via data.table NSE to the subsetting expression and
  ## override the function argument. `head` on a data.table respects the
  ## integer row count directly.
  n_rows <- min(n, nrow(show))
  if (n_rows > 0L) {
    cat("\n")
    print(utils::head(show, n_rows))
  }
  invisible(x)
}
