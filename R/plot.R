#' Plot a `survatr_result`
#'
#' @description
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
#' @family survatr_result methods
#' @examples
#' set.seed(2)
#' n_id <- 40L
#' K <- 5L
#' pp <- data.frame(
#'   id = rep(seq_len(n_id), each = K),
#'   t = rep(seq_len(K), times = n_id),
#'   A = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   Y = rbinom(n_id * K, 1L, 0.1)
#' )
#' fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
#' res <- contrast(
#'   fit,
#'   interventions = list(a1 = causatr::static(1), a0 = causatr::static(0)),
#'   times = 1:5,
#'   type = "survival"
#' )
#' plot(res)
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
    ## Prefer the contrasts view for contrast-kind estimands, but fall back to
    ## curves when the contrasts table is empty -- a single-intervention
    ## quantile is a valid request (a lone median) with no pairwise contrast, so
    ## it must plot the tau-vs-q curve rather than abort.
    which <- if (estimand_is_contrast(x$type) && nrow(x$contrasts) > 0L) {
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
  tbl <- data.table::copy(tbl)
  group_col <- if (identical(which, "contrasts")) "contrast" else "intervention"
  value_col <- if (identical(which, "contrasts")) {
    "estimate"
  } else {
    estimand_field(x$type, "point_col")
  }

  ## Competing-risks results carry a `cause` dimension: draw one line per
  ## (group, cause) pair and label it accordingly. The combined `.plot_group`
  ## column collapses to the plain group label when there is no (non-NA) cause.
  if ("cause" %in% names(tbl) && any(!is.na(tbl[["cause"]]))) {
    tbl[,
      .plot_group := paste0(
        get(group_col),
        " | cause ",
        get("cause")
      )
    ]
    group_col <- ".plot_group"
  }

  groups <- unique(tbl[[group_col]])
  n_g <- length(groups)
  if (is.null(col)) {
    col <- grDevices::palette()[seq_len(min(n_g, length(grDevices::palette())))]
    if (length(col) < n_g) {
      col <- rep_len(col, n_g)
    }
  }

  ## Result index on the x-axis: `time` for the curve estimands, `q` for the
  ## quantile (which plots tau_q against the quantile level).
  idx_col <- estimand_field(x$type, "index")
  if (is.null(xlab)) {
    xlab <- idx_col
  }
  if (is.null(ylab)) {
    ylab <- estimand_field(x$type, "ylab")
  }

  ## y range: point estimates plus CI if present.
  y_vals <- tbl[[value_col]]
  if (ribbon && !all(is.na(tbl$ci_lower))) {
    y_vals <- c(y_vals, tbl$ci_lower, tbl$ci_upper)
  }
  y_range <- range(y_vals, na.rm = TRUE, finite = TRUE)

  plot(
    NA,
    xlim = range(tbl[[idx_col]], na.rm = TRUE),
    ylim = y_range,
    main = main,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  ## Reference line for contrasts: 0 for differences, 1 for ratios (the
  ## null-effect level). Only meaningful on the contrasts view -- a raw
  ## per-intervention curve (e.g. tau_q vs q, a positive time scale) has no
  ## null-effect level, so a contrast-kind estimand plotted as curves gets none.
  op <- estimand_field(x$type, "op")
  if (identical(which, "contrasts") && identical(op, "difference")) {
    graphics::abline(h = 0, lty = 3, col = "grey50")
  } else if (identical(which, "contrasts") && identical(op, "ratio")) {
    graphics::abline(h = 1, lty = 3, col = "grey50")
  }

  for (g_ix in seq_along(groups)) {
    g <- groups[g_ix]
    ## Build the row mask in plain R, then filter, rather than
    ## `tbl[get(group_col) == g]`: evaluated inside [.data.table the symbol
    ## `g` would resolve to a column named "g" if one existed, shadowing the
    ## loop variable. Computing the logical vector outside the subscript can't
    ## be shadowed.
    keep <- tbl[[group_col]] == g
    rows <- tbl[keep]
    data.table::setorderv(rows, idx_col)
    x_vals <- rows[[idx_col]]
    if (ribbon && !all(is.na(rows$ci_lower))) {
      ribbon_col <- grDevices::adjustcolor(col[g_ix], alpha.f = ribbon_alpha)
      ## Closed-polygon CI ribbon: trace the lower bound left-to-right, then
      ## the upper bound right-to-left (`rev()`). The reversal joins the two
      ## bounds into one closed ring so `polygon()` fills the band between
      ## them; concatenating them in the same direction would instead draw a
      ## self-crossing bowtie.
      graphics::polygon(
        x = c(x_vals, rev(x_vals)),
        y = c(rows$ci_lower, rev(rows$ci_upper)),
        col = ribbon_col,
        border = NA
      )
    }
    graphics::lines(
      x_vals,
      rows[[value_col]],
      col = col[g_ix],
      lwd = 2
    )
    graphics::points(x_vals, rows[[value_col]], col = col[g_ix], pch = 19)
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
