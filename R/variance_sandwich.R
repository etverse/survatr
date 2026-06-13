#' Fill sandwich variance into a `survatr_result`
#'
#' @description
#' Take the chunk-2 `estimates` / `contrasts` tables plus the
#' per-intervention IF matrices, aggregate them into pointwise SEs and Wald
#' CIs at `conf_level`, and replace the `NA_real_` placeholder columns.
#'
#' @details
#' Pointwise variance comes from the cross-time IF covariance
#' `V = crossprod(IF) / n_ids^2` (per-individual IFs are stacked rows, so the
#' empirical-mean variance carries the `1 / n_ids^2` scaling). SEs are the
#' square root of its diagonal; Wald CIs are `point +/- z * se` with
#' `z = qnorm(1 - (1 - conf_level) / 2)`. The six contrast types map onto
#' this as:
#'
#' - `survival` / `risk`: SE on `s_hat` / `1 - s_hat` is the same
#'   (`IF_risk = -IF_S`).
#' - `rmst`: `SE(RMST(t_j)) = sqrt(w_j' V w_j)` where `V = crossprod(IF_S) /
#'   n_ids^2` and `w_j` is the trapezoidal-weight row from `rmst_weights()`.
#' - `risk_difference`: contrast IF = `IF_ref - IF_a1`; cross-product / n^2.
#' - `rmst_difference`: contrast IF = `(IF_ref - IF_a1) %*% t(W)`; diagonal
#'   of its cross-product is the pointwise RMST-difference variance.
#' - `risk_ratio`: log-RR IF built pointwise per time; CI computed on the
#'   log scale and exponentiated so the reported `ci_lower` / `ci_upper`
#'   are always strictly positive.
#'
#' Source: M-estimation sandwich aggregation of the cumulative-product
#' survival IF, grounded in Hernán & Robins (2020), *Causal Inference: What
#' If*, Ch. 17; the cross-time delta derivation lives in
#' `CHUNK_3_SANDWICH_A.md`.
#'
#' @param estimates Per-intervention estimates `data.table` from chunk 2.
#' @param contrasts Contrast `data.table` from chunk 2.
#' @param if_list Named list of `compute_survival_if_matrix()` outputs, one
#'   per intervention (names match `estimates$intervention`).
#' @param type Contrast type.
#' @param reference Reference intervention name, or `NULL`.
#' @param times The user time grid.
#' @param conf_level Scalar in (0, 1).
#' @param n_ids Number of individuals.
#'
#' @return A list `list(estimates, contrasts)` with SE / CI columns filled.
#' @noRd
fill_sandwich_ses <- function(
  estimates,
  contrasts,
  if_list,
  type,
  reference,
  times,
  conf_level,
  n_ids
) {
  estimates <- data.table::copy(estimates)
  contrasts <- data.table::copy(contrasts)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  ## --- per-intervention SE (for the `estimates` table) ------------------
  for (iv_name in names(if_list)) {
    IF_mat <- if_list[[iv_name]]$IF_mat
    ## The IF matrix must be n_ids x |times| for `crossprod(IF) / n_ids^2`
    ## to be the cross-time covariance. `compute_survival_if_matrix()`
    ## guarantees this, but assert it so a future regression there surfaces
    ## here rather than as a silently wrong vcov.
    if (
      !is.matrix(IF_mat) ||
        nrow(IF_mat) != n_ids ||
        ncol(IF_mat) != length(times)
    ) {
      rlang::abort(
        paste0(
          "influence-function matrix for intervention \"",
          iv_name,
          "\" is ",
          nrow(IF_mat),
          " x ",
          ncol(IF_mat),
          " but ",
          n_ids,
          " x ",
          length(times),
          " (n_ids x |times|) was expected."
        ),
        class = "survatr_if_failed"
      )
    }
    ## The per-intervention SE family + point column come from the estimand
    ## registry (`estimand_registry.R`). `"pointwise"` reads the diagonal of the
    ## time-covariance directly; `"rmst"` first maps the IF through the
    ## trapezoidal weight matrix (RMST and RMTL share it -- the constant
    ## restriction time has zero gradient, so Var(RMTL) = Var(RMST)).
    level_se <- estimand_field(type, "level_se")
    point_col <- estimand_field(type, "point_col")
    if (identical(level_se, "rmst")) {
      W <- rmst_weights(times)
      IF_se <- IF_mat %*% t(W) ## n_ids x |t|
    } else {
      IF_se <- IF_mat
    }
    se_vec <- sqrt(pmax(diag(crossprod(IF_se)) / n_ids^2, 0))
    estimates[get("intervention") == iv_name, se := se_vec]
    pt_vec <- estimates[get("intervention") == iv_name, get(point_col)]
    estimates[
      get("intervention") == iv_name,
      `:=`(
        ci_lower = pt_vec - z * se_vec,
        ci_upper = pt_vec + z * se_vec
      )
    ]
  }

  ## --- pairwise contrasts (for the `contrasts` table) -------------------
  if (estimand_is_curve(type) || nrow(contrasts) == 0L) {
    return(list(estimates = estimates, contrasts = contrasts))
  }

  other_names <- setdiff(names(if_list), reference)
  ref_S_if <- if_list[[reference]]$IF_mat
  ref_s_hat <- if_list[[reference]]$s_hat
  ## Contrast SE family from the registry: "difference" (plain IF difference),
  ## "logratio" (delta on the log scale, exponentiated CI), or "rmst" (the
  ## trapezoidal quadratic form, shared by RMST and RMTL differences).
  contrast_se <- estimand_field(type, "contrast_se")

  for (a1_name in other_names) {
    a1_S_if <- if_list[[a1_name]]$IF_mat
    a1_s_hat <- if_list[[a1_name]]$s_hat
    ref_risk <- 1 - ref_s_hat
    a1_risk <- 1 - a1_s_hat

    if (identical(contrast_se, "difference")) {
      ## IF on (risk_a1 - risk_a0) = -(IF_S_a1 - IF_S_a0) = IF_S_a0 - IF_S_a1.
      IF_diff <- ref_S_if - a1_S_if
      se_vec <- sqrt(pmax(diag(crossprod(IF_diff)) / n_ids^2, 0))
      est_vec <- contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        estimate
      ]
      contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        `:=`(
          se = se_vec,
          ci_lower = est_vec - z * se_vec,
          ci_upper = est_vec + z * se_vec
        )
      ]
    } else if (identical(contrast_se, "logratio")) {
      ## Delta method on the log scale:
      ##   log(RR(t)) = log(risk_a1(t)) - log(risk_a0(t))
      ##   IF_{log RR}(t) = (1/risk_a1(t)) * IF_risk_a1 - (1/risk_a0(t)) * IF_risk_a0
      ##                  = -(1/risk_a1(t)) * IF_S_a1 + (1/risk_a0(t)) * IF_S_a0
      ## We build IF_log_rr as an n_ids x |t| matrix column-by-column since
      ## the scaling is time-dependent.
      IF_log_rr <- sweep(ref_S_if, 2L, ref_risk, "/") -
        sweep(a1_S_if, 2L, a1_risk, "/")
      ## Guard against division-by-zero when the reference risk is exactly
      ## 0 at some time -- emit NA SE rather than Inf.
      if (any(ref_risk == 0) || any(a1_risk == 0)) {
        bad <- which(ref_risk == 0 | a1_risk == 0)
        IF_log_rr[, bad] <- NA_real_
      }
      se_log <- sqrt(pmax(diag(crossprod(IF_log_rr)) / n_ids^2, 0))
      rr_vec <- contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        estimate
      ]
      ## Report the `se` column on the NATURAL (risk-ratio) scale via the
      ## delta method (multiply the log-scale SE by RR), so it means "SE of
      ## the reported estimand" and matches the bootstrap path (which reports
      ## the SD of the RR replicates). The CI is still built on the log scale
      ## and exponentiated -- transform-respecting and strictly positive -- so
      ## it is intentionally not a symmetric `point +/- z * se` interval.
      se_nat <- rr_vec * se_log
      contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        `:=`(
          se = se_nat,
          ci_lower = exp(log(rr_vec) - z * se_log),
          ci_upper = exp(log(rr_vec) + z * se_log)
        )
      ]
    } else if (identical(contrast_se, "rmst")) {
      ## RMST / RMTL difference: the trapezoidal quadratic form. RMTL difference
      ## = -(RMST difference); the variance is identical (the sign cancels in the
      ## quadratic form), so the same IF and weight matrix serve both and the CI
      ## is anchored on the already-computed `estimate` column.
      W <- rmst_weights(times)
      IF_diff_S <- ref_S_if - a1_S_if
      IF_diff_RMST <- IF_diff_S %*% t(W) ## n_ids x |t|
      se_vec <- sqrt(pmax(diag(crossprod(IF_diff_RMST)) / n_ids^2, 0))
      est_vec <- contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        estimate
      ]
      contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        `:=`(
          se = se_vec,
          ci_lower = est_vec - z * se_vec,
          ci_upper = est_vec + z * se_vec
        )
      ]
    }
  }

  list(estimates = estimates, contrasts = contrasts)
}
