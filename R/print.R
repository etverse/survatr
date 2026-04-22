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
