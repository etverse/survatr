#' Prepare the IPW treatment-model correction (intervention-independent pieces)
#'
#' @description
#' Build the pieces of the stacked-estimating-equation sandwich that come from
#' the **treatment model** under `estimator = "ipw"`: the cross-derivative
#' `A_beta_alpha` of the weighted hazard-MSM score with respect to the
#' propensity coefficients, and the treatment-model IF prep (bread + score).
#' These do not depend on the intervention or the time grid, so are computed
#' once and reused across all interventions, mirroring how
#' `prepare_sandwich_shared()` shares the hazard-model prep.
#'
#' @details
#' The IPW survival estimator is a two-stage M-estimator: stage 1 fits the
#' propensity `alpha` on the baseline rows (one per id); stage 2 fits the
#' weighted hazard MSM `beta` on the at-risk person-period rows with weights
#' `w_i(alpha)`. Because stage-1 `alpha` does not depend on stage-2 `beta`, the
#' stacked bread is block-lower-triangular and the per-id IF on `beta` is
#'
#' ```
#' IF_beta_i = B_bb^{-1} ( psi_beta_i - A_ba B_aa^{-1} psi_alpha_i )
#' ```
#'
#' The first term is the chunk-3 "hazard-only" IF; the second is the
#' treatment-model correction, **subtracted**. This routine builds `A_ba`
#' (= `A_beta_alpha`) by numerically differentiating the mean hazard-MSM score
#' with respect to `alpha` -- exactly the `phi_bar` / `numDeriv::jacobian`
#' construction causatr uses in its point-IPW sandwich
#' (`compute_ipw_if_self_contained_one()`). Only `w_i(alpha)` varies with
#' `alpha`; `mu_eta`, the working residual, and the variance are at `beta_hat`.
#'
#' This block holds because the stabilized weights are **intervention-free** (a
#' single per-id weight at the observed treatment). A future chunk whose weights
#' depend on a quantity that itself depends on `beta` would make the
#' off-diagonal `A_ab` non-zero and break the lower-triangularity.
#'
#' @param fit A `survatr_fit` with `estimator == "ipw"` (carries
#'   `treatment_model`, `marginal_model`, `trim_threshold`).
#' @param prep The hazard-model `causatr:::prepare_model_if()` output (so the
#'   MSM design `X_fit`, sized to the at-risk PP rows, is reused).
#' @param fit_idx Integer `which(fit_rows)` into the sorted `fit$pp_data`.
#' @param id_vec Vector of ids, one per PP row, in `fit$pp_data` order.
#' @param pp_work The working copy of `fit$pp_data` (post `build_risk_set()`),
#'   used to reconstruct the baseline rows.
#'
#' @returns A list with `A_beta_alpha` (`p_beta x p_alpha`) and `prep_trt` (the
#'   treatment-model `prepare_model_if()` output at `n_total = n_ids`).
#' @family variance
#' @noRd
prepare_ipw_correction <- function(fit, prep, fit_idx, id_vec, pp_work) {
  full_tm <- fit$treatment_model
  id <- fit$id
  treatment <- fit$treatment

  ## Reconstruct the baseline rows (one per id) the treatment model was fit on.
  ## `pp_work` is the same data (sorted by id, time) the fit used, so the
  ## first-row-per-id selection reproduces the exact baseline frame.
  baseline_idx <- pp_work[, .I[1L], by = c(id)]$V1
  baseline <- pp_work[baseline_idx]
  baseline_fit <- baseline[full_tm$fit_rows]
  a_obs <- baseline_fit[[treatment]]
  baseline_fit_ids <- baseline_fit[[id]]

  ## Fixed stabilization numerator f(A_i) at the observed treatment, captured
  ## once. Omitting its estimation from the IF is conservative (Robins 1999).
  fixed_numerator <- causatr:::evaluate_density(
    fit$marginal_model,
    a_obs,
    baseline_fit
  )

  weight_closure <- make_survival_weight_closure(
    X_prop = full_tm$X_prop,
    a_obs = a_obs,
    fixed_numerator = fixed_numerator,
    trim_threshold = fit$trim_threshold
  )

  ## Beta-side pieces at beta_hat on the PP fit rows. `prep$X_fit` and the
  ## fitted model's linear predictors / response share the at-risk row order.
  msm <- fit$model
  fam <- msm$family
  eta <- as.numeric(msm$linear.predictors)
  mu <- fam$linkinv(eta)
  mu_eta <- fam$mu.eta(eta)
  var_mu <- fam$variance(mu)
  y_fit <- stats::model.response(stats::model.frame(msm))
  r_fit <- as.numeric(y_fit) - mu
  X_msm <- prep$X_fit
  n_fit <- nrow(X_msm)

  ## Map each PP fit row to its id's position in the baseline fit frame, so the
  ## per-id weight broadcasts onto the hazard score correctly.
  id_map_fit <- match(id_vec[fit_idx], baseline_fit_ids)
  if (anyNA(id_map_fit)) {
    rlang::abort(
      c(
        "IPW sandwich: a person-period fit row has no baseline weight row.",
        i = "The treatment-model id set must match the person-period id set."
      ),
      class = "survatr_ipw_id_mismatch"
    )
  }

  ## Mean hazard-MSM score as a function of alpha (only the weight varies).
  ## phi_bar(alpha) = (1/n_fit) X_msm' diag(w(alpha) mu_eta r / V) 1.
  phi_bar <- function(alpha) {
    w_id <- weight_closure(alpha)
    w_fit <- w_id[id_map_fit]
    s_row <- w_fit * mu_eta * r_fit / var_mu
    as.numeric(crossprod(X_msm, s_row)) / n_fit
  }

  ## Negative-Hessian convention: A_{beta,alpha} = -(1/n) sum d psi / d alpha,
  ## and numDeriv::jacobian(phi_bar) = +(1/n) sum d psi / d alpha, so flip sign.
  ## alpha_hat is the aliased-dropped coefficient vector, aligned with X_prop.
  alpha_hat <- full_tm$alpha_hat
  A_beta_alpha <- -numDeriv::jacobian(phi_bar, x = alpha_hat)

  ## Treatment-model bread + score, scaled by the number of UNITS (= ids), not
  ## the PP row count. `n_total = n_ids` is what makes causatr's
  ## `apply_model_correction()` emit the correction at the same n_ids scaling as
  ## the chunk-3 hazard block, so `crossprod(IF) / n_ids^2` is consistent.
  n_ids <- length(full_tm$fit_rows)
  prep_trt <- causatr:::prepare_model_if(
    model = full_tm$model,
    fit_idx = which(full_tm$fit_rows),
    n_total = n_ids
  )

  list(A_beta_alpha = A_beta_alpha, prep_trt = prep_trt)
}

#' Treatment-model correction matrix for the IPW survival sandwich
#'
#' @description
#' Given the intervention-specific population derivative `J_bar_mat` of the
#' survival curve with respect to the hazard-MSM coefficients, return the
#' `n_ids x |times|` correction matrix to **add** to the chunk-3 hazard-only IF
#' matrix. The correction carries the minus sign of the stacked-EE formula and a
#' single factor of `n_ids` (supplied by `apply_model_correction()`'s
#' `n_total = n_ids`).
#'
#' @details
#' For each time `t`:
#'
#' ```
#' h_t      = n_fit * B_bb^{-1} J_bar(t)          (= A_bb^{-1} J_bar(t))
#' g_alpha  = A_beta_alpha' h_t                   (p_alpha gradient)
#' corr_t   = apply_model_correction(prep_trt, g_alpha)$correction  (n_ids long)
#' TC[, t]  = - corr_t
#' ```
#'
#' The `n_fit` in `h_t` cancels the `1/n_fit` baked into `A_beta_alpha`'s
#' `phi_bar`, leaving the single `n_ids` from `apply_model_correction()`. This
#' mirrors causatr's `compute_ipw_if_self_contained_one()`, with the survival
#' cross-time `J_bar(t)` replacing causatr's scalar marginal-mean gradient.
#'
#' @param J_bar_mat `p_beta x |times|` matrix of population derivatives of
#'   `S^a(t)` wrt the hazard-MSM coefficients (built in
#'   `compute_survival_if_matrix()`).
#' @param B_inv The hazard-model bread `(X'WX)^{-1}` (`prep$B_inv`).
#' @param n_fit Number of at-risk PP rows the hazard MSM was fit on
#'   (`nrow(prep$X_fit)`).
#' @param ipw_corr The `prepare_ipw_correction()` output
#'   (`A_beta_alpha`, `prep_trt`).
#' @param n_ids Number of individuals (rows of the returned matrix).
#'
#' @returns An `n_ids x |times|` numeric matrix to add to the hazard-only IF.
#' @family variance
#' @noRd
compute_treatment_correction <- function(
  J_bar_mat,
  B_inv,
  n_fit,
  ipw_corr,
  n_ids
) {
  A_beta_alpha <- ipw_corr$A_beta_alpha
  prep_trt <- ipw_corr$prep_trt
  n_t <- ncol(J_bar_mat)

  TC <- matrix(0, n_ids, n_t)
  for (j in seq_len(n_t)) {
    ## A_bb^{-1} J_bar(t): the M-estimation bread-projected survival gradient.
    h_t <- n_fit * as.numeric(B_inv %*% J_bar_mat[, j])
    ## g_alpha = A_{beta,alpha}' h_t: the propensity-direction sensitivity.
    g_alpha <- as.numeric(crossprod(A_beta_alpha, h_t))
    ## Per-id propensity correction; subtracted per the block-lower-triangular
    ## M-estimation solution.
    prop_res <- causatr:::apply_model_correction(prep_trt, g_alpha)
    TC[, j] <- -prop_res$correction
  }
  TC
}
