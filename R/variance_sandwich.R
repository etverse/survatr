#' Fill sandwich variance into a `survatr_result`
#'
#' Takes the chunk-2 `estimates` / `contrasts` tables and the list of
#' per-intervention IF matrices, computes pointwise SEs via
#' `crossprod(IF) / n_ids^2`, Wald CIs at `conf_level`, and replaces the
#' `NA_real_` columns. Handles the six contrast types:
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
    vcov_mat <- crossprod(IF_mat) / n_ids^2
    if (type %in% c("survival", "risk", "risk_difference", "risk_ratio")) {
      se_vec <- sqrt(pmax(diag(vcov_mat), 0))
      target_col <- if (type == "survival") "s_hat" else "risk_hat"
      estimates[get("intervention") == iv_name, se := se_vec]
      point_vec <- estimates[get("intervention") == iv_name, get(target_col)]
      estimates[
        get("intervention") == iv_name,
        `:=`(
          ci_lower = point_vec - z * se_vec,
          ci_upper = point_vec + z * se_vec
        )
      ]
    } else if (type %in% c("rmst", "rmst_difference")) {
      W <- rmst_weights(times)
      IF_rmst <- IF_mat %*% t(W) ## n_ids x |t|
      vcov_rmst <- crossprod(IF_rmst) / n_ids^2
      se_vec <- sqrt(pmax(diag(vcov_rmst), 0))
      estimates[get("intervention") == iv_name, se := se_vec]
      rmst_vec <- estimates[get("intervention") == iv_name, get("rmst_hat")]
      estimates[
        get("intervention") == iv_name,
        `:=`(
          ci_lower = rmst_vec - z * se_vec,
          ci_upper = rmst_vec + z * se_vec
        )
      ]
    }
  }

  ## --- pairwise contrasts (for the `contrasts` table) -------------------
  if (type %in% c("survival", "risk", "rmst") || nrow(contrasts) == 0L) {
    return(list(estimates = estimates, contrasts = contrasts))
  }

  other_names <- setdiff(names(if_list), reference)
  ref_S_if <- if_list[[reference]]$IF_mat
  ref_s_hat <- if_list[[reference]]$s_hat

  for (a1_name in other_names) {
    a1_S_if <- if_list[[a1_name]]$IF_mat
    a1_s_hat <- if_list[[a1_name]]$s_hat
    ref_risk <- 1 - ref_s_hat
    a1_risk <- 1 - a1_s_hat

    if (type == "risk_difference") {
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
    } else if (type == "risk_ratio") {
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
      contrasts[
        get("contrast") == paste0(a1_name, " vs ", reference),
        `:=`(
          se = se_log,
          ci_lower = exp(log(rr_vec) - z * se_log),
          ci_upper = exp(log(rr_vec) + z * se_log)
        )
      ]
    } else if (type == "rmst_difference") {
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
