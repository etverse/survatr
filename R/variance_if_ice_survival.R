#' Build a minimal `causatr_fit` for the Track B ICE sandwich
#'
#' @description
#' Hand-construct the smallest `causatr_fit` (type `"longitudinal"`) that the
#' Track B influence-function assembler needs, **without** calling
#' `causatr::fit_ice()` (which would resolve a pseudo family and run causatr's
#' plain, non-survival `ice_iterate()`). The IF machinery only reads a small
#' set of fields off `fit` (`$data`, `$id`, `$time`, `$treatment`,
#' `$details$time_points` / `$n_times` / `$weights` / `$stratified`), so a
#' minimal hand-built object is sufficient and avoids dragging in causatr's
#' contrast path.
#'
#' @param fit A `survatr_fit` with `track == "B"`.
#' @param details Output of `build_ice_surv_details()`.
#' @param data_lag Person-period data with lag + `.survatr_prev_*` columns
#'   (from `prepare_track_b_base()`); becomes `fit$data` for the chain.
#'
#' @return A `causatr_fit` object (type `"longitudinal"`, `model = NULL`).
#' @noRd
build_min_causatr_fit_b <- function(fit, details, data_lag) {
  causatr:::new_causatr_fit(
    model = NULL,
    data = data_lag,
    treatment = fit$treatment,
    outcome = fit$outcome,
    confounders = fit$confounders,
    confounders_tv = fit$confounders_tv,
    family = details$family_outcome,
    estimator = "gcomp",
    type = "longitudinal",
    estimand = "mean",
    id = fit$id,
    time = fit$time,
    censoring = fit$censoring,
    history = fit$history,
    numerator = NULL,
    weights_obj = NULL,
    match_obj = NULL,
    call = fit$call,
    ## Only these `details` fields are read by `ice_if_setup()` and the
    ## survatr survival chain: the time grid, its length, external weights,
    ## and the (NULL = pooled) stratification column.
    details = list(
      time_points = details$time_points,
      n_times = details$n_times,
      weights = details$weights,
      stratified = NULL
    )
  )
}

#' Per-period all-cause event indicators keyed by id
#'
#' @description
#' Build, for each period in the time grid, a named vector mapping individual
#' id (as character) to the all-cause event indicator `D_k` observed at that
#' period. Used by the survival IF chain's failure carry-forward factor
#' `(1 - D_k)` (see `survatr_ice_surv_chain()`).
#'
#' @param data Person-period `data.table` (rectangular; one row per (id, t)).
#' @param time_points Sorted time grid.
#' @param id_col,time_col,outcome Column names.
#'
#' @return A list indexed by step `k`; element `[[k]]` is a named numeric of
#'   `D` at `time_points[k]`, names = ids.
#' @noRd
build_event_by_step <- function(data, time_points, id_col, time_col, outcome) {
  lapply(time_points, function(t_k) {
    rows_k <- data[[time_col]] == t_k
    stats::setNames(
      as.numeric(data[[outcome]][rows_k]),
      as.character(data[[id_col]][rows_k])
    )
  })
}

#' Survival-aware Channel-2 nuisance-correction chain
#'
#' @description
#' The survatr replacement for `causatr:::variance_if_ice_chain()` in the
#' Track B (survival) setting. It mirrors causatr's forward sensitivity
#' cascade -- per step it builds the cross-step gradient `g_k`, calls
#' `causatr:::correct_model()` to get the model-`k` correction and the updated
#' per-individual sensitivity `d_k` -- but injects the survival failure
#' carry-forward factor that causatr's chain (built for a single terminal
#' outcome) lacks.
#'
#' @details
#' The survival-tail pseudo-outcome is `Ỹ_k = D_k + (1 - D_k) q_{k+1}`, so the
#' cross-step derivative is `dỸ_k/dβ_{k+1} = (1 - D_k) dq_{k+1}/dβ_{k+1}` --
#' an extra `(1 - D_k)` versus causatr's plain `dq_{k+1}/dβ_{k+1}`. That factor
#' is **not** absorbed by the at-risk fit-set restriction: at step `k` the fit
#' set (`fit_ids[[k]]`) is at-risk-at-`k`, which **includes** individuals who
#' have the event at `k` (`D_k = 1`). causatr's chain sums the sensitivity over
#' that set unweighted, spuriously propagating the next model's uncertainty
#' through individuals who already failed -- inflating the sandwich, and
#' increasingly so at later horizons (more accumulated failures). This was
#' caught empirically: causatr's verbatim chain over-covered relative to the
#' bootstrap, with the gap growing in `t`.
#'
#' The fix is to multiply each previous-step contribution by `(1 - D_{k})` via
#' `event_by_step`. Everything else -- the bread, the score (`r_score` read off
#' each fitted GLM), the block-triangular back-substitution -- is causatr's,
#' reused through `correct_model()`. Stochastic interventions are rejected for
#' Track B v1, so the Monte-Carlo branches of causatr's chain are omitted here.
#'
#' @param ctx Context list from `causatr:::ice_if_setup()`.
#' @param models_by_step Per-step fitted GLMs (from the `ice_result`).
#' @param fit_ids_by_step Per-step at-risk fit-id character vectors.
#' @param event_by_step Per-step id-keyed `D_k` from `build_event_by_step()`.
#'
#' @return Numeric vector of length `ctx$n`, the accumulated Channel-2 IF.
#' @noRd
survatr_ice_surv_chain <- function(
  ctx,
  models_by_step,
  fit_ids_by_step,
  event_by_step
) {
  data_iv <- ctx$data_iv
  time_points <- ctx$time_points
  n_times <- ctx$n_times
  id_col <- ctx$id_col
  time_col <- ctx$time_col
  all_ids <- ctx$all_ids
  n <- ctx$n
  id_to_idx <- ctx$id_to_idx
  target <- ctx$target
  w_t <- ctx$w_t
  sum_w_target <- ctx$sum_w_target
  has_weights <- ctx$has_weights
  w_at_step <- ctx$w_at_step

  IF_acc <- numeric(n)
  d_vec <- rep(0, n)

  for (step_i in seq_len(n_times)) {
    model_k <- models_by_step[[step_i]]
    if (is.null(model_k)) {
      next
    }
    fit_ids_k <- fit_ids_by_step[[step_i]]
    if (length(fit_ids_k) == 0L) {
      next
    }

    current_time <- time_points[step_i]
    rows_iv_current <- data_iv[[time_col]] == current_time
    iv_data_current <- data_iv[rows_iv_current]
    iv_ids_current <- as.character(iv_data_current[[id_col]])

    X_star_k <- causatr:::iv_design_matrix(model_k, iv_data_current)
    eta_star <- as.numeric(X_star_k %*% causatr:::coef_clean(model_k))
    mu_eta_star <- model_k$family$mu.eta(eta_star)

    if (step_i == 1L) {
      ## g_1 = (1/sum_w) X*' diag(w_target) mu'(eta*): gradient of the ICE mean
      ## w.r.t. the earliest model's parameters (the standardisation step;
      ## everyone is at risk at the first period, no failure factor needed).
      target_in_iv <- match(all_ids[target], iv_ids_current)
      valid_target <- !is.na(target_in_iv)
      target_in_iv <- target_in_iv[valid_target]
      target_w <- w_t[valid_target]
      g_k <- as.numeric(
        crossprod(
          X_star_k[target_in_iv, , drop = FALSE],
          target_w * mu_eta_star[target_in_iv]
        )
      ) /
        sum_w_target
    } else {
      ## g_k = sum_{j in prev fit} w_{k-1,j} (1 - D_{k-1,j}) d_{k-1,j} X*_{k,j}
      ##       mu'_{k,j}. The (1 - D_{k-1,j}) factor is the survival fix.
      prev_fit_ids <- fit_ids_by_step[[step_i - 1L]]
      idx_in_all <- id_to_idx[prev_fit_ids]
      rows_in_iv <- match(prev_fit_ids, iv_ids_current)
      keep <- !is.na(idx_in_all) & !is.na(rows_in_iv)
      if (any(keep)) {
        d_prev <- d_vec[idx_in_all[keep]]
        w_prev <- if (has_weights) {
          unname(w_at_step[[step_i - 1L]][prev_fit_ids[keep]])
        } else {
          rep(1, sum(keep))
        }
        ## Failure carry-forward: zero out individuals who had the event at the
        ## previous step (they do not propagate the next model's sensitivity).
        d_prev_event <- as.numeric(
          event_by_step[[step_i - 1L]][prev_fit_ids[keep]]
        )
        surv_factor <- 1 - d_prev_event
        weights_g <- w_prev *
          surv_factor *
          d_prev *
          mu_eta_star[rows_in_iv[keep]]
        g_k <- as.numeric(
          crossprod(X_star_k[rows_in_iv[keep], , drop = FALSE], weights_g)
        )
      } else {
        g_k <- rep(0, length(causatr:::coef_clean(model_k)))
      }
    }

    fit_id_idx <- id_to_idx[fit_ids_k]
    na_act_k <- model_k$na.action
    if (!is.null(na_act_k)) {
      fit_id_idx <- fit_id_idx[-na_act_k]
    }
    res <- causatr:::correct_model(model_k, g_k, fit_id_idx, n)
    IF_acc <- IF_acc + res$correction
    d_vec <- res$d
  }

  IF_acc
}

#' Per-individual IF for one Track B intervention / horizon
#'
#' @description
#' Channel 1 (sampling term on the baseline survival-tail pseudo-outcome) via
#' `causatr:::ice_if_setup()`, plus Channel 2 (the survival-aware nuisance
#' chain) via `survatr_ice_surv_chain()`. Returns the per-individual influence
#' function for the cumulative risk `R^d(τ)` (the horizon baked into
#' `ice_result$pseudo_final`).
#'
#' @param min_fit Minimal `causatr_fit` from `build_min_causatr_fit_b()`.
#' @param ice_result One horizon's `ice_result` (from
#'   `run_ice_survival_horizon()`).
#' @param target Logical vector over first-period individuals.
#' @param event_by_step Per-step id-keyed `D_k` from `build_event_by_step()`.
#'
#' @return Numeric length-`n` IF on `R^d(τ)`.
#' @noRd
survatr_ice_surv_if_one <- function(
  min_fit,
  ice_result,
  target,
  event_by_step
) {
  ctx <- causatr:::ice_if_setup(min_fit, ice_result, target)
  ctx$IF_vec +
    survatr_ice_surv_chain(
      ctx,
      ice_result$models,
      ice_result$fit_ids,
      event_by_step
    )
}

#' Per-individual IF matrix for the Track B survival curve
#'
#' @description
#' Bridge the per-horizon `ice_result`s (from `compute_ice_survival_curve()`)
#' to the survival sandwich. For each horizon `t`, assemble the per-individual
#' influence function for the cumulative risk `R^d(t)` via
#' `survatr_ice_surv_if_one()`, then negate it to the survival scale
#' (`IF_S = -IF_R`). Stack the columns into the `n_ids x |times|` matrix the
#' shared `fill_sandwich_ses()` consumes.
#'
#' The result is validated numerically against the empirical bootstrap (and an
#' independent delicatessen M-estimation oracle) in `test-ice-survival*.R`.
#'
#' @param min_fit Minimal `causatr_fit` from `build_min_causatr_fit_b()`.
#' @param ice_results Per-horizon list of `ice_result`s.
#' @param times User horizon grid (column count / order).
#' @param target Logical vector over first-period individuals (Track B v1: all
#'   `TRUE`).
#' @param event_by_step Per-step id-keyed `D_k` from `build_event_by_step()`.
#'
#' @return A list with `s_hat` (length `|times|`, survival curve) and `IF_mat`
#'   (`n_ids x |times|` per-individual IF on `S^d(t)`), matching the
#'   `compute_survival_if_matrix()` contract `fill_sandwich_ses()` expects.
#' @noRd
compute_ice_survival_if_matrix <- function(
  min_fit,
  ice_results,
  times,
  target,
  event_by_step
) {
  n <- length(target)
  n_t <- length(times)
  IF_mat <- matrix(0, n, n_t)
  s_hat <- numeric(n_t)

  for (j in seq_len(n_t)) {
    res <- ice_results[[j]]
    if_risk <- survatr_ice_surv_if_one(min_fit, res, target, event_by_step)
    ## Survival scale: S = 1 - R, so the IF flips sign. `na.rm` excludes
    ## entry-censored ids (NA pseudo) from the standardisation mean, matching
    ## the `target` restriction in Channel 1 (Issue #1, 2026-06-03 review).
    IF_mat[, j] <- -if_risk
    s_hat[j] <- 1 - mean(res$pseudo_final, na.rm = TRUE)
  }

  list(s_hat = s_hat, IF_mat = IF_mat)
}
