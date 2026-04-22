#' Forest-plot generic
#'
#' S3 generic for forest-style plots at a user-chosen reference time.
#' survatr ships `forrest.survatr_result()` which slices the result at a
#' single `t_ref` and renders one row per pairwise contrast with a point
#' estimate and CI.
#'
#' The deliberate misspelling `forrest` matches causatr's convention
#' (clarity about the plot style, not the Gump reference).
#'
#' @param x Object to dispatch on.
#' @param ... Arguments passed to methods.
#'
#' @return Method-dependent.
#' @export
forrest <- function(x, ...) {
  UseMethod("forrest")
}

#' Forest plot of contrasts at a reference time
#'
#' Extract the row(s) of `x$contrasts` at `t_ref` and render a
#' horizontal forest plot: one row per contrast, with point estimate,
#' 95% CI (from `x$ci_method`), and a textual label. Only available for
#' contrast-shaped `type`s (`risk_difference`, `risk_ratio`,
#' `rmst_difference`).
#'
#' @param x A `survatr_result`.
#' @param t_ref Numeric scalar; must be in `x$time_grid`.
#' @param col Single color for point + segment. Default
#'   `grDevices::palette()[1]`.
#' @param main Optional plot title.
#' @param xlab Optional x-axis label.
#' @param ... Passed to `plot.default()`.
#'
#' @return The `survatr_result`, invisibly.
#' @method forrest survatr_result
#' @export
forrest.survatr_result <- function(
  x,
  t_ref,
  col = NULL,
  main = NULL,
  xlab = NULL,
  ...
) {
  if (!x$type %in% c("risk_difference", "risk_ratio", "rmst_difference")) {
    rlang::abort(
      paste0(
        "`forrest()` requires a contrast-shaped result (",
        "`type` in {risk_difference, risk_ratio, rmst_difference}). ",
        "Got type = \"",
        x$type,
        "\"."
      ),
      class = "survatr_forrest_wrong_type"
    )
  }
  if (
    !is.numeric(t_ref) ||
      length(t_ref) != 1L ||
      is.na(t_ref) ||
      !(t_ref %in% x$time_grid)
  ) {
    rlang::abort(
      paste0(
        "`t_ref` must be a single value in `x$time_grid` (",
        paste(x$time_grid, collapse = ", "),
        "). Got ",
        deparse(t_ref),
        "."
      ),
      class = "survatr_bad_t_ref"
    )
  }
  rows <- x$contrasts[get("time") == t_ref]
  if (nrow(rows) == 0L) {
    rlang::abort(
      "No contrast rows found at `t_ref`.",
      class = "survatr_bad_t_ref"
    )
  }

  if (is.null(col)) {
    col <- grDevices::palette()[1L]
  }
  if (is.null(main)) {
    main <- paste0("Forest: contrasts at t = ", t_ref)
  }
  if (is.null(xlab)) {
    xlab <- switch(
      x$type,
      risk_difference = "risk difference",
      risk_ratio = "risk ratio",
      rmst_difference = "RMST difference"
    )
  }

  n_c <- nrow(rows)
  y_pos <- seq_len(n_c)
  has_ci <- !all(is.na(rows$ci_lower))
  x_range <- if (has_ci) {
    range(c(rows$estimate, rows$ci_lower, rows$ci_upper), na.rm = TRUE)
  } else {
    range(rows$estimate, na.rm = TRUE)
  }
  ## Add a small margin so points are not clipped on the plot edges.
  x_pad <- 0.05 * diff(x_range)
  x_range <- x_range + c(-x_pad, x_pad)

  graphics::par(mar = c(4.5, 10, 3, 2))
  plot(
    NA,
    xlim = x_range,
    ylim = c(0.5, n_c + 0.5),
    yaxt = "n",
    main = main,
    xlab = xlab,
    ylab = "",
    ...
  )
  ref_line <- if (identical(x$type, "risk_ratio")) 1 else 0
  graphics::abline(v = ref_line, lty = 3, col = "grey50")
  graphics::axis(
    2,
    at = y_pos,
    labels = rows$contrast,
    las = 1,
    cex.axis = 0.9
  )
  if (has_ci) {
    graphics::segments(
      x0 = rows$ci_lower,
      x1 = rows$ci_upper,
      y0 = y_pos,
      y1 = y_pos,
      col = col,
      lwd = 2
    )
  }
  graphics::points(rows$estimate, y_pos, col = col, pch = 19, cex = 1.3)
  invisible(x)
}
