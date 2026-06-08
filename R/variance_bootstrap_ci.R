#' Fill bootstrap SE / CI columns into an existing result
#'
#' Uses the replicate matrix `t` from `bootstrap_survival()` to compute
#' per-column SEs (sample SD) and CIs (percentile or Wald). Updates the
#' `estimates` and `contrasts` data.tables in `res` in place-by-copy.
#'
#' The column layout is intervention-major, cause-mid, time-minor (cause-mid
#' collapses away for single-event types, where `boot$meta$causes` is `NULL`).
#'
#' @param estimates,contrasts data.tables from the point-estimate path.
#' @param boot List from `bootstrap_survival()`.
#' @param type Contrast type.
#' @param conf_level Confidence level.
#' @param boot_ci Either `"percentile"` or `"wald"`.
#'
#' @returns List `list(estimates, contrasts)`.
#' @noRd
fill_bootstrap_ses <- function(
  estimates,
  contrasts,
  boot,
  type,
  conf_level,
  boot_ci
) {
  estimates <- data.table::copy(estimates)
  contrasts <- data.table::copy(contrasts)

  meta <- boot$meta
  alpha <- 1 - conf_level
  z <- stats::qnorm(1 - alpha / 2)

  ## Per-column SE / CI from the replicate matrix.
  col_se <- apply(boot$t, 2L, stats::sd, na.rm = TRUE)
  col_ci <- switch(
    boot_ci,
    percentile = apply(boot$t, 2L, function(v) {
      stats::quantile(v, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    }),
    wald = rbind(
      lower = boot$t0 - z * col_se,
      upper = boot$t0 + z * col_se
    )
  )

  estimand_col <- boot_estimand_col(type)
  block <- meta$n_cause * meta$k_t

  ## Key by (key, cause, time) where present so a (key_val, cause) selector
  ## returns rows in time order, aligning with the bootstrap column blocks.
  key_estimates <- intersect(
    c("intervention", "cause", "time"),
    names(estimates)
  )
  data.table::setkeyv(estimates, key_estimates)
  if (nrow(contrasts) > 0L) {
    key_contrasts <- intersect(c("contrast", "cause", "time"), names(contrasts))
    data.table::setkeyv(contrasts, key_contrasts)
  }

  ## Write one (key_val) block's SE / CI into `tbl`, cause slot by cause slot.
  ## `base_offset` is the column offset of this block in `boot$t`.
  fill_block <- function(tbl, key_col, key_val, base_offset, point_col) {
    for (ci in seq_len(meta$n_cause)) {
      col_idx <- (base_offset + (ci - 1L) * meta$k_t + 1L):(base_offset +
        ci * meta$k_t)
      se_vec <- col_se[col_idx]
      ci_low <- col_ci[1L, col_idx]
      ci_high <- col_ci[2L, col_idx]
      sel <- if (is.null(meta$causes)) {
        tbl[[key_col]] == key_val
      } else {
        tbl[[key_col]] == key_val &
          !is.na(tbl[["cause"]]) &
          tbl[["cause"]] == meta$causes[ci]
      }
      tbl[sel, `:=`(se = se_vec, ci_lower = ci_low, ci_upper = ci_high)]
      ## Wald CI is anchored at the original-fit point (preserved in `tbl`)
      ## rather than the resampled `t0`; percentile already uses the replicate
      ## quantiles directly.
      if (identical(boot_ci, "wald")) {
        pt <- tbl[sel, get(point_col)]
        tbl[sel, `:=`(ci_lower = pt - z * se_vec, ci_upper = pt + z * se_vec)]
      }
    }
    invisible(NULL)
  }

  for (j in seq_along(meta$intervention_names)) {
    iv <- meta$intervention_names[j]
    fill_block(estimates, "intervention", iv, (j - 1L) * block, estimand_col)
  }

  if (meta$n_ctr == 0L) {
    return(list(estimates = estimates, contrasts = contrasts))
  }

  ctr_base <- meta$n_iv * block
  for (j in seq_along(meta$contrast_names)) {
    cn <- meta$contrast_names[j]
    off <- ctr_base + (j - 1L) * block
    fill_block(contrasts, "contrast", cn, off, "estimate")
  }

  list(estimates = estimates, contrasts = contrasts)
}
