#' Fill bootstrap SE / CI columns into an existing result
#'
#' Uses the replicate matrix `t` from `bootstrap_survival()` to compute
#' per-column SEs (sample SD) and CIs (percentile or Wald). Updates the
#' `estimates` and `contrasts` data.tables in `res` in place-by-copy.
#'
#' @param estimates,contrasts data.tables from chunk 2.
#' @param boot List from `bootstrap_survival()`.
#' @param type Contrast type.
#' @param conf_level Confidence level.
#' @param boot_ci Either `"percentile"` or `"wald"`.
#'
#' @return List `list(estimates, contrasts)`.
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

  ## Write into the estimates table. Ordering matches flatten_boot_result:
  ## intervention-major, time-minor.
  estimand_col <- switch(
    type,
    survival = "s_hat",
    risk = "risk_hat",
    risk_difference = "risk_hat",
    risk_ratio = "risk_hat",
    rmst = "rmst_hat",
    rmst_difference = "rmst_hat"
  )
  for (j in seq_along(meta$intervention_names)) {
    iv <- meta$intervention_names[j]
    idx <- ((j - 1L) * meta$k_t + 1L):(j * meta$k_t)
    se_vec <- col_se[idx]
    ci_low <- col_ci[1L, idx]
    ci_high <- col_ci[2L, idx]
    point_vec <- estimates[get("intervention") == iv, get(estimand_col)]
    estimates[
      get("intervention") == iv,
      `:=`(
        se = se_vec,
        ci_lower = ci_low,
        ci_upper = ci_high
      )
    ]
    ## For Wald CIs we want CI around the observed point estimate, which
    ## is already what `t0 +/- z * se` gives. For percentile we want the
    ## replicate quantile, also already what `col_ci` gives. But if the
    ## point estimate moved slightly under resampling (bias), we preserve
    ## the original-fit point in `estimates` and anchor the Wald CI
    ## there rather than at `t0`. Adjust only when boot_ci == "wald".
    if (identical(boot_ci, "wald")) {
      estimates[
        get("intervention") == iv,
        `:=`(
          ci_lower = point_vec - z * se_vec,
          ci_upper = point_vec + z * se_vec
        )
      ]
    }
  }

  if (meta$n_ctr == 0L) {
    return(list(estimates = estimates, contrasts = contrasts))
  }

  for (j in seq_along(meta$contrast_names)) {
    cn <- meta$contrast_names[j]
    col_idx <- meta$ctr_cols[
      ((j - 1L) * meta$k_t + 1L):(j * meta$k_t)
    ]
    se_vec <- col_se[col_idx]
    ci_low <- col_ci[1L, col_idx]
    ci_high <- col_ci[2L, col_idx]
    est_vec <- contrasts[get("contrast") == cn, estimate]
    contrasts[
      get("contrast") == cn,
      `:=`(se = se_vec, ci_lower = ci_low, ci_upper = ci_high)
    ]
    if (identical(boot_ci, "wald")) {
      contrasts[
        get("contrast") == cn,
        `:=`(
          ci_lower = est_vec - z * se_vec,
          ci_upper = est_vec + z * se_vec
        )
      ]
    }
  }

  list(estimates = estimates, contrasts = contrasts)
}
