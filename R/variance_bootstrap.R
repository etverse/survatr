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
  seed,
  causes = NULL
) {
  id_col <- fit$id
  ## Build a per-id row-index map once. `split()` keys on as.character().
  id_vals <- fit$pp_data[[id_col]]
  id_to_rows <- split(seq_len(nrow(fit$pp_data)), id_vals)
  unique_ids <- names(id_to_rows)
  n_ids <- length(unique_ids)

  ## Metadata for column layout. The competing-risks CIF estimands add a cause
  ## dimension between intervention and time (intervention-major, cause-mid,
  ## time-minor); `causes = NULL` collapses to the single-event layout, so the
  ## single-event path is byte-for-byte unchanged.
  iv_names <- names(interventions)
  k_t <- length(times)
  n_iv <- length(iv_names)
  n_cause <- if (is.null(causes)) 1L else length(causes)
  has_contrast <- type %in%
    c(
      "risk_difference",
      "risk_ratio",
      "rmst_difference",
      "cif_difference",
      "cif_ratio"
    )
  contrast_names <- if (has_contrast) {
    paste0(setdiff(iv_names, reference), " vs ", reference)
  } else {
    character(0)
  }
  n_ctr <- length(contrast_names)
  block <- n_cause * k_t
  n_cols <- (n_iv + n_ctr) * block

  meta <- list(
    intervention_names = iv_names,
    contrast_names = contrast_names,
    times = times,
    type = type,
    causes = causes,
    n_cause = n_cause,
    est_cols = seq_len(n_iv * block),
    ctr_cols = if (n_ctr > 0L) {
      (n_iv * block + 1L):(n_iv * block + n_ctr * block)
    } else {
      integer(0)
    },
    n_iv = n_iv,
    n_ctr = n_ctr,
    k_t = k_t
  )

  ## Weight handling depends on the estimator. For gcomp with external weights,
  ## the weight vector simply travels with the sampled id's rows (no re-
  ## estimation). For IPW the weights are a FITTED quantity (stabilized
  ## density ratios), so they must be re-estimated inside `surv_fit()` on each
  ## replicate -- that is exactly what makes the bootstrap a full two-stage
  ## resample that propagates the treatment-model uncertainty. Carrying the
  ## original IPW weights would freeze the propensity fit and undercount the
  ## variance, so we do NOT travel them for IPW.
  is_ipw <- identical(fit$estimator, "ipw")
  orig_weights <- if (is_ipw) NULL else fit$weights

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
        estimator = fit$estimator,
        model_fn = fit$model_fn,
        propensity_model_fn = if (is_ipw) {
          fit$propensity_model_fn
        } else {
          stats::glm
        },
        trim = fit$trim,
        ## Track B (ice) needs the time-varying confounders + lag order; both
        ## are NULL for Track A, where `surv_fit()` ignores them.
        confounders_tv = fit$confounders_tv,
        history = if (is.null(fit$history)) Inf else fit$history,
        ## Competing risks: re-activate the cause-specific path each replicate
        ## (NULL for single-event fits). The cause models are re-estimated per
        ## resample, propagating their uncertainty into the bootstrap CI.
        competing = fit$competing,
        ## IPCW: re-estimate the censoring model per replicate. `ipcw = NULL`
        ## when IPCW was not requested; the arguments are ignored by surv_fit()
        ## when NULL. Per-replicate IPCW trim re-quantiles the weights within
        ## the replicate (matching the per-replicate treatment-weight trim).
        ipcw = fit$ipcw,
        censoring_model_fn = if (!is.null(fit$ipcw)) {
          fit$censoring_model_fn
        } else {
          stats::glm
        }
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
        cause = causes,
        ci_method = "none"
      ),
      error = function(e) NULL
    )
    if (is.null(boot_res)) {
      return(rep(NA_real_, n_cols))
    }
    flatten_boot_result(boot_res, meta)
  }

  ## Reproducibility under parallel backends requires an L'Ecuyer
  ## stream: `mclapply` / `parLapply` ignore the serial RNG state
  ## otherwise, so two calls with the same `seed` would produce
  ## different replicate sequences. Save and restore the prior RNG
  ## kind so we do not leak L'Ecuyer back to the caller.
  prev_kind <- NULL
  if (!is.null(seed)) {
    if (!identical(parallel, "no")) {
      prev_kind <- RNGkind()
      RNGkind("L'Ecuyer-CMRG")
    }
    set.seed(seed)
  }
  on.exit(
    if (!is.null(prev_kind)) {
      RNGkind(prev_kind[1L], prev_kind[2L], prev_kind[3L])
    },
    add = TRUE
  )

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

#' Map an estimand column to the bootstrap value column
#'
#' @param type Contrast type.
#'
#' @returns The name of the per-intervention estimand column for `type`.
#' @noRd
boot_estimand_col <- function(type) {
  switch(
    type,
    survival = "s_hat",
    risk = "risk_hat",
    risk_difference = "risk_hat",
    risk_ratio = "risk_hat",
    rmst = "rmst_hat",
    rmst_difference = "rmst_hat",
    rmtl = "rmtl_hat",
    rmtl_difference = "rmtl_hat",
    cif = "cif_hat",
    cif_difference = "cif_hat",
    cif_ratio = "cif_hat"
  )
}

#' Flatten a `survatr_result` to the bootstrap's column vector
#'
#' Column layout, intervention-major then cause-mid then time-minor:
#' - First `n_iv * n_cause * k_t` entries: per-intervention point estimates
#'   (`s_hat` for `survival`, `risk_hat` for `risk` / `risk_difference` /
#'   `risk_ratio`, `rmst_hat` for `rmst*`, `cif_hat` for `cif*`).
#' - Next `n_ctr * n_cause * k_t` entries: contrast estimates, same ordering.
#'
#' For single-event types `n_cause == 1` (`meta$causes` is `NULL`) and the
#' layout collapses to the original intervention-major, time-minor order.
#'
#' @param res A `survatr_result` from `contrast(..., ci_method = "none")`.
#' @param meta The metadata list from the caller's layout.
#'
#' @returns Numeric vector of length `(n_iv + n_ctr) * n_cause * k_t`.
#' @noRd
flatten_boot_result <- function(res, meta) {
  estimand_col <- boot_estimand_col(meta$type)
  block <- meta$n_cause * meta$k_t

  ## Extract one (intervention | contrast) block: cause-major, time-minor.
  ## `meta$causes = NULL` => a single cause slot with no cause filter, i.e. the
  ## single-event layout.
  pull <- function(tbl, key_col, key_val, value_col) {
    out <- numeric(block)
    for (ci in seq_len(meta$n_cause)) {
      rows <- if (is.null(meta$causes)) {
        tbl[get(key_col) == key_val]
      } else {
        tbl[get(key_col) == key_val & get("cause") == meta$causes[ci]]
      }
      data.table::setkeyv(rows, "time")
      idx <- ((ci - 1L) * meta$k_t + 1L):(ci * meta$k_t)
      out[idx] <- rows[[value_col]]
    }
    out
  }

  est_vec <- numeric(meta$n_iv * block)
  for (j in seq_along(meta$intervention_names)) {
    iv <- meta$intervention_names[j]
    idx <- ((j - 1L) * block + 1L):(j * block)
    est_vec[idx] <- pull(res$estimates, "intervention", iv, estimand_col)
  }

  if (meta$n_ctr == 0L) {
    return(est_vec)
  }
  ctr_vec <- numeric(meta$n_ctr * block)
  for (j in seq_along(meta$contrast_names)) {
    cn <- meta$contrast_names[j]
    idx <- ((j - 1L) * block + 1L):(j * block)
    ctr_vec[idx] <- pull(res$contrasts, "contrast", cn, "estimate")
  }
  c(est_vec, ctr_vec)
}
