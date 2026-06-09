#' Print a `survatr_fit`
#'
#' @description
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
#' @family survatr_fit functions
#' @examples
#' set.seed(1)
#' n_id <- 30L
#' K <- 4L
#' pp <- data.frame(
#'   id = rep(seq_len(n_id), each = K),
#'   t = rep(seq_len(K), times = n_id),
#'   A = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   Y = rbinom(n_id * K, 1L, 0.1)
#' )
#' fit <- surv_fit(pp, "Y", "A", ~1, "id", "t", time_formula = ~ factor(t))
#' print(fit)
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
  ## Track B (ICE) carries time-varying confounders + a Markov lag order that
  ## Track A does not; surface them so the printout reflects the actual model.
  if (identical(x$track, "B")) {
    tv <- if (is.null(x$confounders_tv)) {
      "none"
    } else {
      paste(all.vars(x$confounders_tv), collapse = " + ")
    }
    hist_lab <- if (is.null(x$history) || is.infinite(x$history)) {
      "full"
    } else {
      as.character(x$history)
    }
    cat(sprintf("  TV covars:   %s\n", tv))
    cat(sprintf("  History:     %s\n", hist_lab))
  }
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
#' @description
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
#'   type = "risk_difference"
#' )
#' print(res)
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

  ## Competing-risks CIF contrasts condition on surviving the competing events;
  ## repeat the caveat here so a printed result never shows the numbers silently.
  if (x$type %in% c("cif_difference", "cif_ratio")) {
    cat(
      "  Note:        cause-specific CIF contrasts condition on surviving the\n",
      "               competing events (truncation by death); see the\n",
      "               competing-risks vignette.\n",
      sep = ""
    )
  }

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

#' Print a `survatr_diag`
#'
#' @description
#' Compact per-panel summary for the diagnostic object returned by
#' `diagnose()`. Reports how many periods have positivity flags, the range
#' of SMDs in the balance panel, IPW weight summary (if applicable), and the
#' competing-risks identity-check result (if applicable).
#'
#' @param x A `survatr_diag`.
#' @param ... Unused.
#'
#' @returns The diagnostic object, invisibly.
#' @family survatr_fit functions
#' @examples
#' set.seed(4)
#' n_id <- 40L; K <- 3L
#' pp <- data.frame(
#'   id = rep(seq_len(n_id), each = K),
#'   t  = rep(seq_len(K), times = n_id),
#'   A  = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   L  = rep(rnorm(n_id), each = K),
#'   Y  = rbinom(n_id * K, 1L, 0.05)
#' )
#' fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
#' print(diagnose(fit))
#' @export
print.survatr_diag <- function(x, ...) {
  cat("<survatr_diag>\n")

  ## Positivity panel
  pos <- x$positivity
  if (!is.null(pos) && nrow(pos) > 0L) {
    n_flag_low <- sum(pos$flag_low, na.rm = TRUE)
    n_flag_high <- sum(pos$flag_high, na.rm = TRUE)
    cat(sprintf(
      "  Positivity:  %d periods, hazards [%.4f, %.4f]",
      nrow(pos),
      min(pos$h_min, na.rm = TRUE),
      max(pos$h_max, na.rm = TRUE)
    ))
    if (n_flag_low + n_flag_high > 0L) {
      cat(sprintf(
        " [!] %d low, %d high",
        n_flag_low,
        n_flag_high
      ))
    }
    cat("\n")
  }

  ## Balance panel
  bal <- x$balance
  if (!is.null(bal) && nrow(bal) > 0L) {
    smd_range <- range(bal$smd, na.rm = TRUE)
    cat(sprintf(
      "  Balance:     %d variable(s), SMD range [%.3f, %.3f]\n",
      length(unique(bal$variable)),
      smd_range[1L],
      smd_range[2L]
    ))
  } else {
    cat("  Balance:     (no confounders or not applicable)\n")
  }

  ## Weights panel (IPW only)
  wt <- x$weights
  if (!is.null(wt)) {
    cat(sprintf(
      "  Weights:     ESS = %.1f / %d, max = %.3f, top5%% share = %.3f\n",
      wt$ess,
      wt$n_ids,
      wt$max_weight,
      wt$top5_share
    ))
  }

  ## Censoring panel
  cens <- x$censoring
  if (!is.null(cens) && nrow(cens) > 0L) {
    total_cens <- sum(cens$n_censored, na.rm = TRUE)
    total_pp <- sum(cens$n_pp_rows, na.rm = TRUE)
    cat(sprintf(
      "  Censoring:   %d events over %d person-periods (%.1f%%)\n",
      total_cens,
      total_pp,
      100 * total_cens / max(total_pp, 1L)
    ))
  } else if (is.null(x$censoring)) {
    cat("  Censoring:   (no censoring column supplied)\n")
  }

  ## Competing panel
  cr <- x$competing
  if (!is.null(cr)) {
    id_chk <- attr(cr, "identity_check")
    max_dev <- if (!is.null(id_chk)) {
      max(id_chk$abs_dev, na.rm = TRUE)
    } else {
      NA_real_
    }
    cat(sprintf(
      "  Competing:   %d cause(s), max identity deviation = %.2e\n",
      length(unique(cr$cause[!is.na(cr$cause)])),
      max_dev
    ))
    cat(
      "  Note:        CIFs condition on surviving competing events",
      "(truncation by death).\n"
    )
  }

  invisible(x)
}
