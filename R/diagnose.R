#' Diagnose a fitted survatr hazard model
#'
#' @description
#' Returns per-period diagnostic panels for a `survatr_fit`:
#' \describe{
#'   \item{positivity}{Predicted hazard distribution by time; flags periods
#'     with hazards near 0 or 1.}
#'   \item{balance}{Standardized mean differences (or treatment–confounder
#'     correlation for continuous treatment) by time.}
#'   \item{weights}{IPW only: per-id weight summary (ESS, max, top-5% share).
#'     `NULL` for gcomp / ice.}
#'   \item{censoring}{Per-period censoring incidence by arm. `NULL` when no
#'     censoring column was supplied.}
#'   \item{competing}{Competing risks only: per-cause CIF + all-cause survival
#'     + identity check `Σ F^(j)(t) + S(t) = 1`. `NULL` for single-event fits.}
#' }
#'
#' @details
#' All panels operate on the at-risk rows as defined by `build_risk_set()` —
#' the same rows the hazard model was fit on. For Track B (ICE) the hazard
#' model is fit lazily inside `contrast()`, so the positivity panel reports
#' the empirical per-period event rate rather than model-predicted hazards.
#'
#' @param fit A `survatr_fit` returned by `surv_fit()`.
#' @param ... Unused; reserved for future arguments.
#'
#' @returns A `survatr_diag` object: a named list with elements `positivity`,
#'   `balance`, `weights`, `censoring`, `competing`. Each is a `data.table`
#'   indexed by `time` (and `cause` for `competing`), or `NULL` when the panel
#'   does not apply.
#'
#' @family survatr_fit functions
#' @seealso [surv_fit()], [contrast()]
#' @examples
#' set.seed(3)
#' n_id <- 50L; K <- 4L
#' pp <- data.frame(
#'   id   = rep(seq_len(n_id), each = K),
#'   t    = rep(seq_len(K), times = n_id),
#'   A    = rep(rbinom(n_id, 1L, 0.5), each = K),
#'   L    = rep(rnorm(n_id), each = K),
#'   Y    = rbinom(n_id * K, 1L, 0.05)
#' )
#' fit <- surv_fit(pp, "Y", "A", ~L, "id", "t", time_formula = ~1)
#' dx  <- diagnose(fit)
#' print(dx)
#' @export
diagnose <- function(fit, ...) UseMethod("diagnose")

#' @rdname diagnose
#' @export
diagnose.survatr_fit <- function(fit, ...) {
  pp <- data.table::copy(fit$pp_data)

  pos_dt <- diag_positivity(fit, pp)
  bal_dt <- diag_balance(fit, pp)
  wt_dt <- if (identical(fit$estimator, "ipw")) diag_weights(fit) else NULL
  cens_dt <- if (!is.null(fit$censoring)) diag_censoring(fit, pp) else NULL
  cr_dt <- if (!is.null(fit$competing)) diag_competing(fit, pp) else NULL

  new_survatr_diag(
    positivity = pos_dt,
    balance = bal_dt,
    weights = wt_dt,
    censoring = cens_dt,
    competing = cr_dt
  )
}

#' Constructor for `survatr_diag`
#'
#' @param positivity,balance,weights,censoring,competing Panel data.tables (or
#'   `NULL`). See [diagnose.survatr_fit()] for schema.
#'
#' @returns A list of class `survatr_diag`.
#' @noRd
new_survatr_diag <- function(
  positivity,
  balance,
  weights,
  censoring,
  competing
) {
  structure(
    list(
      positivity = positivity,
      balance = balance,
      weights = weights,
      censoring = censoring,
      competing = competing
    ),
    class = "survatr_diag"
  )
}

## ── Panel helpers ────────────────────────────────────────────────────────────

#' Per-period hazard positivity panel
#'
#' Predicts the hazard on the at-risk rows and summarises by period. For Track B
#' (ICE), the per-step models are fit lazily inside `contrast()`; the empirical
#' event rate among at-risk individuals is used instead.
#'
#' @param fit A `survatr_fit`.
#' @param pp Copy of `fit$pp_data`.
#'
#' @returns `data.table` with columns `time | n_at_risk | h_min | h_mean |
#'   h_max | flag_low | flag_high`.
#' @noRd
diag_positivity <- function(fit, pp) {
  fit_rows <- build_risk_set(pp, fit$outcome, fit$id, fit$censoring)
  pp_ar <- pp[fit_rows]

  if (!is.null(fit$cause_models)) {
    ## Competing-risks: all-cause hazard from summed cause-specific predictions.
    ## `cr_augment_pp` adds `.cf_all_haz` column.
    pp_aug <- data.table::copy(pp_ar)
    pp_aug <- cr_augment_pp(
      pp_aug,
      fit$cause_models,
      fit$causes,
      fit$id,
      fit$time
    )
    h_obs <- pp_aug[[".cf_all_haz"]]
  } else if (!is.null(fit$model)) {
    ## gcomp / ipw: model-predicted hazard on observed covariates.
    h_obs <- as.numeric(
      stats::predict(fit$model, newdata = pp_ar, type = "response")
    )
  } else {
    ## ICE (Track B): no fitted model yet; use empirical event rate.
    h_obs <- as.numeric(pp_ar[[fit$outcome]])
  }

  pp_ar[, .survatr_h := h_obs]
  pos_dt <- pp_ar[,
    .(
      n_at_risk = .N,
      h_min = min(.survatr_h, na.rm = TRUE),
      h_mean = mean(.survatr_h, na.rm = TRUE),
      h_max = max(.survatr_h, na.rm = TRUE),
      ## Flags fire when any individual-period hazard exceeds the positivity
      ## bounds -- a signal that the model may be extrapolating in sparse cells.
      flag_low = any(.survatr_h < 0.001, na.rm = TRUE),
      flag_high = any(.survatr_h > 0.999, na.rm = TRUE)
    ),
    by = c(fit$time)
  ]
  data.table::setnames(pos_dt, fit$time, "time")
  pos_dt[]
}

#' Per-period covariate balance panel
#'
#' Computes the standardized mean difference (SMD) for binary treatment, or the
#' treatment–confounder correlation for continuous treatment, at each period
#' among at-risk individuals. Both `confounders` (baseline) and
#' `confounders_tv` (Track B time-varying) are included.
#'
#' SMD = (mean(L|A=1) - mean(L|A=0)) / pooled SD. For continuous treatment,
#' a Pearson correlation `cor(A, L)` per time is returned in the `smd` column.
#'
#' @param fit A `survatr_fit`.
#' @param pp Copy of `fit$pp_data`.
#'
#' @returns `data.table` with columns `time | variable | smd | n_a1 | n_a0`
#'   (n_a0 = NA for continuous treatment).
#' @noRd
diag_balance <- function(fit, pp) {
  base_vars <- all.vars(fit$confounders)
  tv_vars <- if (!is.null(fit$confounders_tv)) {
    all.vars(fit$confounders_tv)
  } else {
    character(0L)
  }
  all_vars <- union(base_vars, tv_vars)
  all_vars <- intersect(all_vars, names(pp)) # keep only present columns

  if (length(all_vars) == 0L) {
    return(data.table::data.table())
  }

  fit_rows <- build_risk_set(pp, fit$outcome, fit$id, fit$censoring)
  pp_ar <- pp[fit_rows]
  trt_col <- fit$treatment

  ## Determine whether to compute SMD (binary) or correlation (continuous).
  trt_vals <- unique(pp_ar[[trt_col]])
  is_binary <- length(na.omit(trt_vals)) == 2L &&
    all(na.omit(trt_vals) %in% c(0L, 1L, 0, 1))

  rows <- lapply(fit$time_grid, function(tk) {
    pp_t <- pp_ar[get(fit$time) == tk]
    lapply(all_vars, function(v) {
      if (is_binary) {
        vals1 <- pp_t[get(trt_col) == 1L][[v]]
        vals0 <- pp_t[get(trt_col) == 0L][[v]]
        m1 <- mean(vals1, na.rm = TRUE)
        m0 <- mean(vals0, na.rm = TRUE)
        s1 <- stats::sd(vals1, na.rm = TRUE)
        s0 <- stats::sd(vals0, na.rm = TRUE)
        ps <- sqrt((s1^2 + s0^2) / 2)
        smd <- if (!is.na(ps) && ps > 0) (m1 - m0) / ps else NA_real_
        data.table::data.table(
          time = tk,
          variable = v,
          smd = smd,
          n_a1 = length(vals1),
          n_a0 = length(vals0)
        )
      } else {
        ## Pearson correlation as a continuous-treatment balance measure.
        smd <- tryCatch(
          stats::cor(pp_t[[trt_col]], pp_t[[v]], use = "complete.obs"),
          error = function(e) NA_real_
        )
        data.table::data.table(
          time = tk,
          variable = v,
          smd = smd,
          n_a1 = NA_integer_,
          n_a0 = NA_integer_
        )
      }
    })
  })
  data.table::rbindlist(unlist(rows, recursive = FALSE))[]
}

#' IPW per-id weight distribution panel
#'
#' Extracts the stabilized density-ratio weights from the fit (one per id,
#' broadcast across periods in `fit$weights`). Reports the effective sample
#' size (ESS = (Σw)² / Σw²), maximum weight, and share of total weight held
#' by the top 5% of units.
#'
#' @param fit A `survatr_fit` with `estimator = "ipw"`.
#'
#' @returns `data.table` with one row: `ess | max_weight | top5_share | n_ids`.
#' @noRd
diag_weights <- function(fit) {
  pp <- fit$pp_data
  ## Under Track A IPW the stabilized weight is constant within id (baseline
  ## weight broadcast to all PP rows). Extract one value per id via the first
  ## row of each id in the PP data.
  first_idx <- pp[, .I[1L], by = c(fit$id)][[2L]]
  w_ids <- fit$weights[first_idx]

  sw <- sum(w_ids)
  sw2 <- sum(w_ids^2)
  ess <- if (sw2 > 0) sw^2 / sw2 else NA_real_
  top5 <- w_ids >= stats::quantile(w_ids, 0.95, names = FALSE)
  top5_share <- if (sw > 0) sum(w_ids[top5]) / sw else NA_real_

  data.table::data.table(
    ess = ess,
    max_weight = max(w_ids),
    top5_share = top5_share,
    n_ids = length(w_ids)
  )
}

#' Per-period censoring incidence panel
#'
#' Computes, for each period and each treatment arm, the number and proportion
#' of at-risk individuals who were censored. A large arm imbalance in the
#' censoring rate signals potential informative censoring that motivates IPCW.
#'
#' @param fit A `survatr_fit` with a non-`NULL` `censoring` column.
#' @param pp Copy of `fit$pp_data`.
#'
#' @returns `data.table` with columns `time | arm | n_at_risk | n_censored |
#'   prop_censored`.
#' @noRd
diag_censoring <- function(fit, pp) {
  cens_col <- fit$censoring
  trt_col <- fit$treatment
  time_col <- fit$time

  ## Use the full PP grid (before risk-set filtering) to see all censoring
  ## events, including the one at the period when censoring occurs.
  cens_dt <- pp[,
    .(
      n_at_risk = .N,
      n_censored = sum(get(cens_col) == 1L, na.rm = TRUE),
      prop_censored = mean(get(cens_col) == 1L, na.rm = TRUE)
    ),
    by = c(time_col, trt_col)
  ]
  data.table::setnames(cens_dt, c(time_col, trt_col), c("time", "arm"))
  cens_dt[]
}

#' Competing-risks CIF panel
#'
#' Computes the marginal (observed-treatment) per-cause CIF and all-cause
#' survival at each time point, then checks the partition-of-unity identity
#' `Σ_j F^(j)(t) + S(t) = 1`. Reports the maximum absolute deviation from 1.
#'
#' The caveat that cause-specific CIFs condition on surviving the competing
#' events (truncation by death) is attached as the `caveat` attribute.
#'
#' @param fit A competing-risks `survatr_fit` (non-`NULL` `cause_models`).
#' @param pp Copy of `fit$pp_data`.
#'
#' @returns `data.table` with columns `cause | time | cif_hat` plus a
#'   `data.table` attribute `identity_check` (one row per time: `time |
#'   sum_cif_plus_surv | abs_dev`).
#' @noRd
diag_competing <- function(fit, pp) {
  ## Augment with per-cause and all-cause hazards using the marginal
  ## (observed-treatment) PP data.
  pp_aug <- cr_augment_pp(pp, fit$cause_models, fit$causes, fit$id, fit$time)

  times <- fit$time_grid
  cif_dt <- compute_cif_curves(
    pp_aug,
    fit$id,
    fit$time,
    times,
    fit$causes,
    "observed"
  )
  s_dt <- compute_cr_survival_curve(pp_aug, fit$id, fit$time, times, "observed")

  ## Keep only the columns needed for the panel.
  cif_panel <- cif_dt[, .(cause, time, cif_hat)]

  ## Identity check: at each time, the sum of all CIFs + all-cause survival
  ## should equal 1 (up to floating-point). A large deviation indicates a
  ## numerical issue in the curve computation.
  cif_sum <- cif_dt[, .(sum_cif = sum(cif_hat)), by = "time"]
  surv_dt <- s_dt[, .(time, s_hat)]
  id_check <- cif_sum[surv_dt, on = "time"]
  id_check[, sum_cif_plus_surv := sum_cif + s_hat]
  id_check[, abs_dev := abs(sum_cif_plus_surv - 1)]
  id_check[, c("sum_cif", "s_hat") := NULL]

  data.table::setattr(
    cif_panel,
    "identity_check",
    id_check[]
  )
  ## Truncation-by-death caveat: cause-specific CIFs condition on surviving
  ## the competing causes; the CIF captures cause-j probability among those
  ## who have not yet experienced any competing event.
  data.table::setattr(
    cif_panel,
    "caveat",
    paste0(
      "Cause-specific CIFs condition on surviving competing events ",
      "(truncation by death). The CIF for cause j is the probability of ",
      "cause-j failure among those not yet failed from any cause."
    )
  )
  cif_panel[]
}
