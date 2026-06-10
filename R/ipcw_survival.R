## IPCW (per-period censoring weights) for the IPW weighted hazard MSM.
## Separated from `ipw_survival.R` to keep each file under the ~300-line cap.
## causatr supplies the IF primitives; survatr builds the per-period running
## product and the stacked-EE censoring-model block.

#' Fit the censoring hazard model on the person-period grid
#'
#' @description
#' Fits the per-period censoring hazard
#' `logit P(C_k = 1 | A, L, at-risk through k) = gamma(k) + delta_A A + delta_L L`
#' on at-risk, not-yet-censored rows (those where neither the event nor a prior
#' censoring has occurred). This is broader than the hazard-MSM fit set: it
#' includes the censoring-event rows themselves (where `C_k = 1`) because the
#' censoring outcome is observed there. The baseline-hazard `gamma(k)` uses
#' `time_formula` (same interface as the outcome hazard). Treatment is always
#' included; additional covariates come from `ipcw_formula`.
#'
#' @param data Person-period `data.table`, already validated and sorted by
#'   `(id, time)`. Must have `.survatr_prev_event` and `.survatr_prev_cens`
#'   columns from a prior `build_risk_set()` call.
#' @param censoring Column name of the censoring indicator (`1` = censored, `0`/`NA` = at risk).
#' @param treatment Column name of the treatment variable.
#' @param ipcw_formula One-sided formula for the censoring-model covariates
#'   (e.g. `~ L1 + L2`). Treatment and the time basis are added automatically;
#'   this formula supplies only the additional covariate adjustment.
#' @param time_formula One-sided time-basis formula (same as `surv_fit()`'s
#'   `time_formula`; the same basis is used for the censoring hazard's
#'   `gamma(k)`).
#' @param censoring_model_fn Fitting function (`stats::glm`-compatible).
#' @param id,time Column names.
#' @param ... Forwarded to `censoring_model_fn` and the numerator model fit
#'   (excluding `family`, which is always `binomial()`).
#'
#' @returns A list with:
#'   \describe{
#'     \item{`cens_model`}{Fitted denominator censoring hazard model.}
#'     \item{`num_model`}{Fitted numerator (treatment + time only) model.}
#'     \item{`cens_fit_rows`}{Logical mask selecting the censoring-model fit rows.}
#'   }
#' @family survatr_fit functions
#' @noRd
fit_censoring_model <- function(
  data,
  censoring,
  treatment,
  ipcw_formula,
  time_formula,
  censoring_model_fn,
  id,
  time,
  ...
) {
  ## Censoring model fit rows: at-risk and not-yet-censored, including the
  ## actual censoring rows (C_k = 1). Differs from the hazard-MSM fit rows
  ## which additionally require `is_uncensored()` (excluding the C_k = 1
  ## rows themselves). The `.survatr_prev_*` columns are written by
  ## `build_risk_set()` before this call.
  cens_fit_rows <- data[[".survatr_prev_event"]] == 0 &
    data[[".survatr_prev_cens"]] == 0
  cens_fit_data <- data[cens_fit_rows]

  ## Build the full denominator formula: C ~ time_basis + A + ipcw_covariates.
  ## `as.character(formula)[2]` extracts the RHS as a string so we can
  ## concatenate with the treatment name and the ipcw covariates without
  ## re-parsing the formula object.
  time_rhs <- as.character(time_formula)[2L]
  ipcw_rhs <- if (identical(ipcw_formula, ~1)) {
    NULL
  } else {
    as.character(ipcw_formula)[2L]
  }
  cens_rhs_parts <- c(time_rhs, treatment, ipcw_rhs)
  cens_rhs <- paste(cens_rhs_parts, collapse = " + ")
  cens_form <- stats::as.formula(paste(censoring, "~", cens_rhs))

  ## Denominator: full model conditioning on treatment, time, and confounders.
  ## `family = binomial()` is forced (censoring is a 0/1 indicator).
  cens_model <- censoring_model_fn(
    formula = cens_form,
    data = cens_fit_data,
    family = stats::binomial(),
    ...
  )

  ## Numerator: treatment + time only — the marginal censoring probability
  ## that stabilizes the cumulative product (Cole & Hernán 2004, Am J Epidemiol
  ## 160:461-468). Omitting the numerator estimation from the IF is conservative;
  ## the denominator-model correction is what narrows the naive weights-as-known SE.
  num_rhs <- paste(time_rhs, "+", treatment)
  num_form <- stats::as.formula(paste(censoring, "~", num_rhs))
  num_model <- censoring_model_fn(
    formula = num_form,
    data = cens_fit_data,
    family = stats::binomial(),
    ...
  )

  list(
    cens_model = cens_model,
    num_model = num_model,
    cens_fit_rows = cens_fit_rows
  )
}

#' Compute per-period cumulative IPCW weights on the full PP grid
#'
#' @description
#' Forms the stabilized inverse-probability-of-censoring weight for each
#' person-period row:
#'
#' ```
#' W^C_{i,k} = prod_{m <= k} (1 - g_num_{i,m}) / (1 - g_{i,m})
#' ```
#'
#' where `g_{i,m}` is the conditional censoring hazard from `cens_model` and
#' `g_num_{i,m}` is the marginal hazard from `num_model`. With
#' `num_model = NULL` the numerator collapses to `1` (unstabilized weight).
#'
#' Weights are computed on ALL PP rows (the running product at row `(i, k)`
#' uses predictions at periods `m = 1, ..., k`). A per-period winsorization at
#' the `trim`-th quantile may be applied (Cole & Hernán 2008). The resolved
#' per-period cutoffs are returned as a named vector so the sandwich can clip
#' at the same fixed thresholds (re-quantiling inside the numerical-derivative
#' closure would make the weights non-smooth in γ).
#'
#' @param data Person-period `data.table`, sorted by `(id, time)`.
#' @param cens_model Fitted denominator censoring hazard model.
#' @param num_model Fitted numerator model, or `NULL` (unstabilized).
#' @param id,time Column names.
#' @param trim Numeric scalar in `(0, 1]` or `NULL`. When `< 1`, weights
#'   at each period are winsorized at the `trim`-th quantile.
#'
#' @returns A list with:
#'   \describe{
#'     \item{`weights`}{Numeric vector, length `nrow(data)`: per-row `W^C_{i,k}`.}
#'     \item{`trim_thresholds`}{Named numeric vector (by time period) of fixed
#'       winsorization cutoffs, or `NULL` when `trim` is `NULL` / `>= 1`.}
#'   }
#' @family survatr_fit functions
#' @noRd
compute_ipcw_running_weights <- function(
  data,
  cens_model,
  num_model = NULL,
  id,
  time,
  trim = NULL
) {
  ## Denominator: g_{i,m} = P(C_m = 1 | ...) for every PP row. Predict on ALL
  ## rows (not just cens_fit_rows) so the running product is well-defined at
  ## rows that were censored or had the event in prior periods.
  g_den <- as.numeric(
    stats::predict(cens_model, newdata = data, type = "response")
  )

  ## Per-row factor: (1 - g_num) / (1 - g_den) for stabilized,
  ## or 1 / (1 - g_den) for unstabilized.
  if (!is.null(num_model)) {
    g_num <- as.numeric(
      stats::predict(num_model, newdata = data, type = "response")
    )
    per_row_factor <- (1 - g_num) / (1 - g_den)
  } else {
    per_row_factor <- 1 / (1 - g_den)
  }

  ## Running product within id (data is already sorted by (id, time) from
  ## prepare_pp_data()). `stats::ave` with FUN = cumprod is the clearest
  ## approach but cumprod is not directly supported; use data.table for
  ## the within-id running product.
  id_vec_all <- data[[id]]
  W_C <- ipcw_running_cumprod(per_row_factor, id_vec_all)

  ## Per-period trim: winsorize at the `trim`-th quantile within each period.
  ## Late-period cumulative products accumulate the most variance inflation,
  ## so trimming per-period (rather than globally) targets the right tail.
  ## The per-period thresholds are fixed here and reused by the sandwich.
  trim_thresholds <- NULL
  if (!is.null(trim) && trim < 1) {
    time_vec <- data[[time]]
    unique_times <- sort(unique(time_vec))
    trim_thresholds <- vapply(
      unique_times,
      function(tt) {
        stats::quantile(W_C[time_vec == tt], trim, names = FALSE)
      },
      numeric(1L)
    )
    names(trim_thresholds) <- as.character(unique_times)
    for (tt in unique_times) {
      mask <- time_vec == tt
      W_C[mask] <- pmin(W_C[mask], trim_thresholds[as.character(tt)])
    }
  }

  list(weights = W_C, trim_thresholds = trim_thresholds)
}

#' Running cumulative product within id (data.table helper)
#'
#' @description
#' Compute `cumprod(v)` within groups defined by `id_vec`, assuming the data is
#' already sorted so that rows of the same id are contiguous (guaranteed by the
#' `(id, time)` key from `prepare_pp_data()`).
#'
#' @param v Numeric vector of per-row factors.
#' @param id_vec Vector of ids, one per row (same order as `v`).
#'
#' @returns Numeric vector of length `length(v)`.
#' @family survatr_fit functions
#' @noRd
ipcw_running_cumprod <- function(v, id_vec) {
  dt_tmp <- data.table::data.table(.f = v, .id = id_vec)
  dt_tmp[, .W := cumprod(.f), by = ".id"]
  dt_tmp[[".W"]]
}

#' Build the gamma-closure for the stabilized IPCW weight
#'
#' @description
#' Return a closure `gamma -> W^C_{i,k}(gamma)` for all MSM fit rows. Used by
#' the sandwich to numerically differentiate the weighted hazard-MSM mean score
#' with respect to the censoring-model coefficients `gamma`. Only the
#' denominator `1 - g_{i,m}(gamma)` varies with `gamma`; the stabilization
#' numerator `(1 - g_num_{i,m})` is captured fixed.
#'
#' The censoring design matrix on ALL PP rows (`X_cens_all`) is pre-built
#' via `causatr:::iv_design_matrix()` so the closure avoids calling
#' `stats::predict()` (which would require rebuilding the model frame on each
#' numDeriv perturbation and is ~10x slower).
#'
#' @param X_cens_all `n_pp x p_cens` censoring design matrix on all PP rows
#'   (from `causatr:::iv_design_matrix(cens_model, pp_data)`).
#' @param numer_factor_all Length-`n_pp` fixed numerator factor per PP row:
#'   `1 - g_num_{i,m}` for stabilized, or `1` for unstabilized.
#' @param id_vec_all PP-row id vector (length `n_pp`).
#' @param fit_idx_msm Integer indices of MSM fit rows into the PP data.
#' @param trim_thresholds Named numeric vector of per-period winsorization
#'   cutoffs (from the point estimate), or `NULL`.
#' @param time_vec_msm Time values at the MSM fit rows (for trim lookup).
#'
#' @returns A function of `gamma` (length `p_cens`) returning a numeric
#'   vector of IPCW weights at the MSM fit rows.
#' @family variance
#' @noRd
make_ipcw_weight_closure_pp <- function(
  X_cens_all,
  numer_factor_all,
  id_vec_all,
  fit_idx_msm,
  trim_thresholds,
  time_vec_msm
) {
  force(X_cens_all)
  force(numer_factor_all)
  force(id_vec_all)
  force(fit_idx_msm)
  force(trim_thresholds)
  force(time_vec_msm)
  has_trim <- !is.null(trim_thresholds)

  function(gamma) {
    eta_C <- as.numeric(X_cens_all %*% gamma)
    g_C <- stats::plogis(eta_C)
    per_row_factor <- numer_factor_all / (1 - g_C)
    W_C_all <- ipcw_running_cumprod(per_row_factor, id_vec_all)
    W_C_msm <- W_C_all[fit_idx_msm]
    if (has_trim) {
      tt_char <- as.character(time_vec_msm)
      caps <- trim_thresholds[tt_char]
      W_C_msm <- pmin(W_C_msm, caps)
    }
    W_C_msm
  }
}

#' Prepare the IPCW censoring-model correction (intervention-independent pieces)
#'
#' @description
#' Build the pieces of the stacked-EE sandwich that come from the
#' **censoring model**: the cross-derivative `A_beta_gamma` of the weighted
#' hazard-MSM mean score with respect to the censoring-model coefficients, and
#' the per-individual censoring-model IF matrix. These are intervention-
#' independent and computed once, shared across all `contrast()` interventions
#' (mirroring `prepare_ipw_correction()` for the treatment block).
#'
#' The censoring-model IF is computed as:
#' ```
#' IF_gamma_i = n_ids * (sum_{k in cens_fit_i} r_score_k * X_cens_k) * (X_cens' X_cens)^{-1}
#' ```
#' where `r_score_k` is the per-row score contribution from
#' `causatr:::prepare_model_if()`. The sum within id accounts for the multiple
#' per-period rows in the censoring model (unlike the one-row-per-id treatment
#' model). Scaling by `n_ids` gives IF_i of magnitude O(1) so
#' `crossprod(IF) / n_ids^2` is the correct sandwich covariance.
#'
#' @param fit A `survatr_fit` with `ipcw` active (carries `censoring_model`,
#'   `ipcw_numerator_model`, `ipcw_trim_thresholds`, `ipw_treatment_weights_pp`).
#' @param prep The hazard-model `causatr:::prepare_model_if()` output.
#' @param fit_idx Integer `which(fit_rows)` into `pp_work`.
#' @param id_vec Vector of ids, one per PP row, in `pp_work` order.
#' @param unique_ids Vector of unique ids.
#' @param pp_work Working copy of `fit$pp_data`.
#'
#' @returns A list with `A_beta_gamma` (`p_beta x p_gamma`) and
#'   `IF_gamma_per_id` (`n_ids x p_gamma`).
#' @family variance
#' @noRd
prepare_ipcw_correction <- function(
  fit,
  prep,
  fit_idx,
  id_vec,
  unique_ids,
  pp_work
) {
  id <- fit$id
  time_col <- fit$time
  n_ids <- length(unique_ids)

  ## Rebuild censoring-model fit rows: at-risk AND not-yet-censored.
  ## `.survatr_prev_*` were stripped from pp_data; we call build_risk_set()
  ## on a fresh copy of pp_work so the `.survatr_prev_*` columns are added.
  pp_cens <- data.table::copy(pp_work)
  build_risk_set(pp_cens, fit$outcome, id, fit$censoring)
  cens_fit_rows <- pp_cens[[".survatr_prev_event"]] == 0 &
    pp_cens[[".survatr_prev_cens"]] == 0
  fit_idx_cens <- which(cens_fit_rows)
  n_cens_fit <- length(fit_idx_cens)

  ## Censoring-model bread + per-row score, at n_total = n_cens_fit (the
  ## standard GLM scaling over the (i,k) fit rows). We aggregate within id
  ## manually below and rescale by n_ids to get the per-individual IF.
  prep_cens <- causatr:::prepare_model_if(
    model = fit$censoring_model,
    fit_idx = fit_idx_cens,
    n_total = n_cens_fit
  )

  ## Per-id summed censoring-model score (n_ids x p_cens).
  ## `B_inv_cens = n_cens_fit * (X_cens' X_cens)^{-1}` and
  ## `A_gamma^{-1} = n_ids * (X_cens' X_cens)^{-1}`, so:
  ## `IF_gamma_per_id = n_ids * psi_per_id %*% (X_cens' X_cens)^{-1}
  ##                  = (n_ids / n_cens_fit) * psi_per_id %*% B_inv_cens`.
  X_cens_fit <- prep_cens$X_fit
  r_cens_fit <- prep_cens$r_score
  id_cens_fit <- id_vec[fit_idx_cens]
  id_cens_f <- factor(id_cens_fit, levels = unique_ids)
  weighted_cens <- X_cens_fit * r_cens_fit
  psi_cens_per_id_raw <- rowsum(weighted_cens, id_cens_f, reorder = FALSE)
  ## Pad rows for ids with no censoring-model fit rows (edge case: all censored
  ## at period 1); those rows contribute zero IF.
  p_cens <- ncol(X_cens_fit)
  psi_cens_per_id <- matrix(0, n_ids, p_cens)
  rownames(psi_cens_per_id) <- as.character(unique_ids)
  psi_cens_per_id[rownames(psi_cens_per_id_raw), ] <- psi_cens_per_id_raw
  IF_gamma_per_id <- (n_ids / n_cens_fit) * psi_cens_per_id %*% prep_cens$B_inv

  ## Pre-build the censoring design matrix on ALL PP rows for the numDeriv
  ## closure (faster than calling predict() on each perturbation).
  X_cens_all <- causatr:::iv_design_matrix(fit$censoring_model, pp_work)

  ## Fixed numerator factors (not a function of gamma): P(C_m = 0 | A, m).
  if (!is.null(fit$ipcw_numerator_model)) {
    g_num_all <- as.numeric(
      stats::predict(
        fit$ipcw_numerator_model,
        newdata = pp_work,
        type = "response"
      )
    )
    numer_factor_all <- 1 - g_num_all
  } else {
    # nocov start
    ## Unstabilized path: numerator = 1 everywhere. Reachable when the public
    ## API exposes stabilize_ipcw = FALSE (not yet in v1); preserved for
    ## forward-compatibility with the compute_ipcw_running_weights contract.
    numer_factor_all <- rep(1, nrow(pp_work))
    # nocov end
  }

  ## Treatment weights at MSM fit rows (fixed in the gamma perturbation;
  ## only the IPCW part varies with gamma).
  w_ipw_msm <- fit$ipw_treatment_weights_pp[fit_idx]

  ## MSM score components at the hazard-MSM fit rows (fixed at beta_hat).
  msm <- fit$model
  fam <- msm$family
  eta_msm <- as.numeric(msm$linear.predictors)
  mu_msm <- fam$linkinv(eta_msm)
  mu_eta_msm <- fam$mu.eta(eta_msm)
  var_msm <- fam$variance(mu_msm)
  y_fit <- stats::model.response(stats::model.frame(msm))
  r_msm <- as.numeric(y_fit) - mu_msm
  X_msm_fit <- prep$X_fit
  n_fit_msm <- nrow(X_msm_fit)

  ## Build the IPCW weight closure: gamma -> W^C_{i,k}(gamma) at MSM fit rows.
  time_vec_msm <- pp_work[[time_col]][fit_idx]
  ipcw_closure <- make_ipcw_weight_closure_pp(
    X_cens_all = X_cens_all,
    numer_factor_all = numer_factor_all,
    id_vec_all = pp_work[[id]],
    fit_idx_msm = fit_idx,
    trim_thresholds = fit$ipcw_trim_thresholds,
    time_vec_msm = time_vec_msm
  )

  ## Mean hazard-MSM score as a function of gamma (only the IPCW weight varies).
  ## `phi_bar_cens(gamma)` = (1/n_fit) X_msm' diag(combined_w(gamma) * mu_eta * r / V) 1.
  ## The treatment weight w_ipw_msm is fixed; only W^C varies with gamma.
  phi_bar_cens <- function(gamma) {
    W_C_msm <- ipcw_closure(gamma)
    combined_w <- w_ipw_msm * W_C_msm
    s_row <- combined_w * mu_eta_msm * r_msm / var_msm
    as.numeric(crossprod(X_msm_fit, s_row)) / n_fit_msm
  }

  ## A_{beta,gamma} = -(1/n) Σ d psi_beta / d gamma, so flip sign of jacobian.
  gamma_hat <- stats::coef(fit$censoring_model)
  A_beta_gamma <- -numDeriv::jacobian(phi_bar_cens, x = gamma_hat)

  list(A_beta_gamma = A_beta_gamma, IF_gamma_per_id = IF_gamma_per_id)
}

#' Censoring-model correction matrix for the IPCW survival sandwich
#'
#' @description
#' Given the intervention-specific population derivative `J_bar_mat` and the
#' pre-built IPCW correction pieces, return the `n_ids x |times|` correction
#' matrix to **add** to the chunk-3 + IPW-treatment IF matrix. The correction
#' propagates the censoring-model score through the cross-derivative and the
#' cross-time delta chain, exactly mirroring `compute_treatment_correction()`
#' for the treatment block.
#'
#' For each time `t`:
#' ```
#' h_t         = n_fit * B_bb^{-1} J_bar(t)
#' g_gamma     = A_beta_gamma' h_t
#' corr_t      = IF_gamma_per_id %*% g_gamma
#' CC[, t]     = -corr_t
#' ```
#'
#' @param J_bar_mat `p_beta x |times|` matrix of population derivatives
#'   of `S^a(t)` w.r.t. the hazard-MSM coefficients.
#' @param B_inv Hazard-MSM bread (`prep$B_inv`).
#' @param n_fit Number of at-risk PP rows the hazard MSM was fit on.
#' @param ipcw_corr Output of `prepare_ipcw_correction()`.
#' @param n_ids Number of individuals.
#'
#' @returns An `n_ids x |times|` numeric matrix to add to the IF matrix.
#' @family variance
#' @noRd
compute_censoring_correction <- function(
  J_bar_mat,
  B_inv,
  n_fit,
  ipcw_corr,
  n_ids
) {
  A_beta_gamma <- ipcw_corr$A_beta_gamma
  IF_gamma_per_id <- ipcw_corr$IF_gamma_per_id
  n_t <- ncol(J_bar_mat)

  CC <- matrix(0, n_ids, n_t)
  for (j in seq_len(n_t)) {
    ## A_bb^{-1} J_bar(t): the M-estimation bread-projected survival gradient.
    h_t <- n_fit * as.numeric(B_inv %*% J_bar_mat[, j])
    ## g_gamma = A_{beta,gamma}' h_t: censoring-direction sensitivity.
    g_gamma <- as.numeric(crossprod(A_beta_gamma, h_t))
    ## Per-id censoring correction; subtracted per the stacked-EE solution.
    CC[, j] <- -as.numeric(IF_gamma_per_id %*% g_gamma)
  }
  CC
}
