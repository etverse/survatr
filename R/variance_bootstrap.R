#' Empirical bootstrap for Track A survival curves / contrasts
#'
#' Resamples **individuals** (all of each id's PP rows together), refits the
#' hazard model on each resample, recomputes the per-intervention survival
#' curves and the requested contrasts, and returns a `B x K` matrix of
#' replicate statistics where `K = n_interventions * |t| + n_contrasts * |t|`
#' for contrast-shaped types or `n_interventions * |t|` for curve-only types.
#'
#' Resampling by id preserves the within-id cumulative-product dependence
#' structure. Per-row resampling would break it and severely understate
#' variance for longer time horizons.
#'
#' @param fit A `survatr_fit`.
#' @param interventions Named list of `causatr_intervention` objects.
#' @param times User time grid (sorted unique, validated upstream).
#' @param type Contrast type.
#' @param reference Reference intervention name, or `NULL`.
#' @param n_boot Integer >= 1, number of replicates.
#' @param parallel One of `"no"`, `"multicore"`, `"snow"` (forwarded to
#'   `boot::boot()`).
#' @param ncpus Integer >= 1.
#' @param seed Integer or `NULL`; when non-null, `set.seed(seed)` before
#'   the bootstrap loop so the entire replicate sequence is reproducible.
#'
#' @return A list with
#' - `t0` -- numeric vector of the point-estimate quantities in the same
#'   order as `t`'s columns (from the original fit);
#' - `t` -- `B x K` matrix of replicate quantities;
#' - `meta` -- a list describing the column layout (`est_cols`,
#'   `ctr_cols`, `intervention_names`, `contrast_names`) so downstream
#'   code can map columns back to estimates / contrasts rows.
#' @noRd
bootstrap_survival <- function(
  fit,
  interventions,
  times,
  type,
  reference,
  n_boot,
  parallel,
  ncpus,
  seed
) {
  id_col <- fit$id
  ## Build a per-id row-index map once. `split()` keys on as.character().
  id_vals <- fit$pp_data[[id_col]]
  id_to_rows <- split(seq_len(nrow(fit$pp_data)), id_vals)
  unique_ids <- names(id_to_rows)
  n_ids <- length(unique_ids)

  ## Metadata for column layout.
  iv_names <- names(interventions)
  k_t <- length(times)
  n_iv <- length(iv_names)
  has_contrast <- type %in%
    c("risk_difference", "risk_ratio", "rmst_difference")
  contrast_names <- if (has_contrast) {
    paste0(setdiff(iv_names, reference), " vs ", reference)
  } else {
    character(0)
  }
  n_ctr <- length(contrast_names)
  n_cols <- n_iv * k_t + n_ctr * k_t

  meta <- list(
    intervention_names = iv_names,
    contrast_names = contrast_names,
    times = times,
    type = type,
    est_cols = seq_len(n_iv * k_t),
    ctr_cols = if (n_ctr > 0L) {
      (n_iv * k_t + 1L):(n_iv * k_t + n_ctr * k_t)
    } else {
      integer(0)
    },
    n_iv = n_iv,
    n_ctr = n_ctr,
    k_t = k_t
  )

  ## Preserve per-id weights when the original fit was weighted. Chunk 4
  ## does not add weight-aware resampling on top of what the user already
  ## supplied -- the weight vector simply travels with the sampled id's
  ## rows. Chunk 5 (IPW) will add density-ratio weights as a separate path.
  orig_weights <- fit$weights

  ## `boot::boot()` passes (data, indices) to the statistic function. We
  ## do not use `data_arg` directly -- the id vector is closed over as
  ## `unique_ids` -- but the signature must match boot's contract.
  statistic_fn <- function(data_arg, indices) {
    sampled_ids <- unique_ids[indices]
    row_blocks <- vector("list", n_ids)
    weight_blocks <- if (!is.null(orig_weights)) vector("list", n_ids) else NULL
    for (i in seq_along(sampled_ids)) {
      rows <- id_to_rows[[sampled_ids[i]]]
      block <- data.table::copy(fit$pp_data[rows])
      ## Re-id each sampled individual with a bootstrap-local integer so
      ## doubled ids do not collapse back together in `prepare_pp_data()`.
      block[, (id_col) := i]
      row_blocks[[i]] <- block
      if (!is.null(weight_blocks)) {
        weight_blocks[[i]] <- orig_weights[rows]
      }
    }
    boot_pp <- data.table::rbindlist(row_blocks)
    boot_w <- if (!is.null(weight_blocks)) unlist(weight_blocks) else NULL

    boot_fit <- tryCatch(
      surv_fit(
        data = boot_pp,
        outcome = fit$outcome,
        treatment = fit$treatment,
        confounders = fit$confounders,
        id = id_col,
        time = fit$time,
        censoring = fit$censoring,
        time_formula = fit$time_formula,
        weights = boot_w,
        estimator = "gcomp",
        model_fn = fit$model_fn
      ),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) {
      return(rep(NA_real_, n_cols))
    }
    boot_res <- tryCatch(
      contrast(
        boot_fit,
        interventions = interventions,
        times = times,
        type = type,
        reference = reference,
        ci_method = "none"
      ),
      error = function(e) NULL
    )
    if (is.null(boot_res)) {
      return(rep(NA_real_, n_cols))
    }
    flatten_boot_result(boot_res, meta)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  b <- boot::boot(
    data = unique_ids,
    statistic = statistic_fn,
    R = n_boot,
    parallel = parallel,
    ncpus = ncpus
  )

  ## Failure guard. boot::boot() does not abort on per-replicate errors;
  ## our statistic_fn returns NA_real_ vectors when a refit / contrast
  ## fails. If more than 10% of replicates failed, abort loudly rather
  ## than quietly producing wide CIs.
  fail_mask <- apply(b$t, 1L, function(row) all(is.na(row)))
  if (mean(fail_mask) > 0.10) {
    rlang::abort(
      paste0(
        "Bootstrap: ",
        sum(fail_mask),
        " of ",
        n_boot,
        " replicates failed to refit / contrast (>10%). Try reducing the ",
        "covariate set, using a larger n, or switching to ",
        "`ci_method = \"sandwich\"`."
      ),
      class = "survatr_boot_failed"
    )
  }

  list(t0 = b$t0, t = b$t, meta = meta, n_failed = sum(fail_mask))
}

#' Flatten a `survatr_result` to the bootstrap's column vector
#'
#' Column layout:
#' - First `n_iv * k_t` entries: per-intervention point estimates
#'   (`s_hat` for `survival`, `risk_hat` for `risk` / `risk_difference`
#'   / `risk_ratio`, `rmst_hat` for `rmst` / `rmst_difference`).
#'   Ordered intervention-major, time-minor (all times for iv1, then all
#'   times for iv2, ...).
#' - Next `n_ctr * k_t` entries: contrast estimates, same ordering.
#'
#' @param res A `survatr_result` from `contrast(..., ci_method = "none")`.
#' @param meta The metadata list from the caller's layout.
#'
#' @return Numeric vector of length `n_iv * k_t + n_ctr * k_t`.
#' @noRd
flatten_boot_result <- function(res, meta) {
  estimand_col <- switch(
    meta$type,
    survival = "s_hat",
    risk = "risk_hat",
    risk_difference = "risk_hat",
    risk_ratio = "risk_hat",
    rmst = "rmst_hat",
    rmst_difference = "rmst_hat"
  )

  est_vec <- numeric(meta$n_iv * meta$k_t)
  for (j in seq_along(meta$intervention_names)) {
    iv <- meta$intervention_names[j]
    rows <- res$estimates[get("intervention") == iv]
    data.table::setkeyv(rows, "time")
    vals <- rows[[estimand_col]]
    idx <- ((j - 1L) * meta$k_t + 1L):(j * meta$k_t)
    est_vec[idx] <- vals
  }

  if (meta$n_ctr == 0L) {
    return(est_vec)
  }
  ctr_vec <- numeric(meta$n_ctr * meta$k_t)
  for (j in seq_along(meta$contrast_names)) {
    cn <- meta$contrast_names[j]
    rows <- res$contrasts[get("contrast") == cn]
    data.table::setkeyv(rows, "time")
    vals <- rows[["estimate"]]
    idx <- ((j - 1L) * meta$k_t + 1L):(j * meta$k_t)
    ctr_vec[idx] <- vals
  }
  c(est_vec, ctr_vec)
}

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
