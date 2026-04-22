#' Plot a `survatr_result`
#'
#' Base-R graphics rendering of survival / risk / RMST curves or of
#' risk-difference / risk-ratio / RMST-difference contrasts. Dispatches
#' on `x$type`:
#'
#' - `survival` / `risk` / `rmst`: one line per intervention (colored),
#'   CI ribbons when populated.
#' - `risk_difference` / `risk_ratio` / `rmst_difference`: one line per
#'   contrast (colored), CI ribbons when populated, and a reference line
#'   at 0 (for differences) or 1 (for ratios).
#'
#' @param x A `survatr_result`.
#' @param which One of `"auto"`, `"curves"`, `"contrasts"`. `"auto"`
#'   picks `curves` for `survival` / `risk` / `rmst` and `contrasts` for
#'   the three pairwise types.
#' @param main Optional plot title.
#' @param xlab,ylab Optional axis labels. Sensible defaults are picked
#'   based on `type`.
#' @param col Optional vector of colors, length `n_interventions` (or
#'   `n_contrasts`). Defaults to `grDevices::palette()`.
#' @param ribbon Logical; when `TRUE` and CI columns are populated, shade
#'   a polygon between `ci_lower` and `ci_upper`. Default `TRUE`.
#' @param ribbon_alpha Numeric in `(0, 1)`; ribbon transparency (via
#'   `grDevices::adjustcolor`). Default `0.2`.
#' @param ... Passed to the underlying `plot.default()` call.
#'
#' @return The `survatr_result`, invisibly.
#' @export
plot.survatr_result <- function(
  x,
  which = c("auto", "curves", "contrasts"),
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  col = NULL,
  ribbon = TRUE,
  ribbon_alpha = 0.2,
  ...
) {
  which <- match.arg(which)
  if (identical(which, "auto")) {
    which <- if (
      x$type %in% c("risk_difference", "risk_ratio", "rmst_difference")
    ) {
      "contrasts"
    } else {
      "curves"
    }
  }
  if (identical(which, "contrasts") && nrow(x$contrasts) == 0L) {
    rlang::abort(
      paste0(
        "`which = \"contrasts\"` requested but `x$contrasts` is empty ",
        "(type = \"",
        x$type,
        "\" has no pairwise contrasts)."
      ),
      class = "survatr_plot_no_contrasts"
    )
  }

  tbl <- if (identical(which, "contrasts")) x$contrasts else x$estimates
  group_col <- if (identical(which, "contrasts")) "contrast" else "intervention"
  value_col <- if (identical(which, "contrasts")) {
    "estimate"
  } else {
    switch(
      x$type,
      survival = "s_hat",
      risk = "risk_hat",
      rmst = "rmst_hat",
      risk_difference = "risk_hat",
      risk_ratio = "risk_hat",
      rmst_difference = "rmst_hat"
    )
  }

  groups <- unique(tbl[[group_col]])
  n_g <- length(groups)
  if (is.null(col)) {
    col <- grDevices::palette()[seq_len(min(n_g, length(grDevices::palette())))]
    if (length(col) < n_g) {
      col <- rep_len(col, n_g)
    }
  }

  if (is.null(xlab)) {
    xlab <- "time"
  }
  if (is.null(ylab)) {
    ylab <- switch(
      x$type,
      survival = "S(t)",
      risk = "1 - S(t)",
      rmst = "RMST(t)",
      risk_difference = "risk difference",
      risk_ratio = "risk ratio",
      rmst_difference = "RMST difference"
    )
  }

  ## y range: point estimates plus CI if present.
  y_vals <- tbl[[value_col]]
  if (ribbon && !all(is.na(tbl$ci_lower))) {
    y_vals <- c(y_vals, tbl$ci_lower, tbl$ci_upper)
  }
  y_range <- range(y_vals, na.rm = TRUE, finite = TRUE)

  plot(
    NA,
    xlim = range(tbl$time, na.rm = TRUE),
    ylim = y_range,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  if (x$type %in% c("risk_difference", "rmst_difference")) {
    graphics::abline(h = 0, lty = 3, col = "grey50")
  } else if (identical(x$type, "risk_ratio")) {
    graphics::abline(h = 1, lty = 3, col = "grey50")
  }

  for (g_ix in seq_along(groups)) {
    g <- groups[g_ix]
    rows <- tbl[get(group_col) == g]
    data.table::setorder(rows, time)
    if (ribbon && !all(is.na(rows$ci_lower))) {
      ribbon_col <- grDevices::adjustcolor(col[g_ix], alpha.f = ribbon_alpha)
      graphics::polygon(
        x = c(rows$time, rev(rows$time)),
        y = c(rows$ci_lower, rev(rows$ci_upper)),
        col = ribbon_col,
        border = NA
      )
    }
    graphics::lines(
      rows$time,
      rows[[value_col]],
      col = col[g_ix],
      lwd = 2
    )
    graphics::points(rows$time, rows[[value_col]], col = col[g_ix], pch = 19)
  }

  graphics::legend(
    "bottomleft",
    legend = as.character(groups),
    col = col,
    lwd = 2,
    bty = "n"
  )
  invisible(x)
}
