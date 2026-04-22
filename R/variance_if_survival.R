#' Compute the per-individual influence function matrix for S^a(t)
#'
#' Implements the delta-method cross-time IF chain described in
#' `CHUNK_3_SANDWICH_A.md`:
#'
#' ```
#' IF_i(t) = (S_i(t) - S_bar(t))                                  ## Ch1
#'         + J_bar(t) %*% B_inv %*% psi_i                         ## Ch2
#' ```
#'
#' where
#'
#' - `S_i(t) = prod_{k <= t} (1 - h^a_{i,k})` is individual i's counterfactual
#'   survival;
#' - `J_bar(t) = -(1/n_ids) * sum_i [S_i(t) * sum_{k <= t} s_{i,k} * X_{i,k}]`
#'   is the population derivative of `S_bar(t)` wrt beta, with per-row
#'   sensitivity `s_{i,k} := mu_eta(eta^a_{i,k}) / (1 - h^a_{i,k})`;
#' - `psi_i = sum_{(i,k) in fit_rows} r_score_{i,k} * X_{i,k}^fit` is the
#'   per-individual score contribution; and
#' - `B_inv` is the bread (from `causatr:::prepare_model_if()`).
#'
#' For the binomial / quasibinomial logit link this simplifies to
#' `s_{i,k} = h_{i,k}` (the `(1 - h)` in the denominator cancels the
#' `h(1 - h)` in `mu_eta`), but the code uses the general formula so a
#' future chunk that swaps in a non-logit link does not need to revisit the
#' derivation.
#'
#' @param fit A `survatr_fit`.
#' @param intervention A single `causatr_intervention`.
#' @param times User-supplied time grid (sorted unique, validated upstream).
#' @param prep The output of `causatr:::prepare_model_if(fit$model, fit_idx,
#'   n_total)` -- shared across interventions for efficiency.
#' @param fit_idx Integer vector: `which(fit_rows)` relative to the sorted
#'   `fit$pp_data`.
#' @param id_vec Character / integer vector of ids, one per PP row, in the
#'   same order as `fit$pp_data`.
#' @param unique_ids Vector of unique id values (in the order they first
#'   appear in `fit$pp_data`).
#'
#' @return A list with `s_hat` (length `|times|`) and `IF_mat`
#'   (`n_ids x |times|` matrix of per-individual IFs on `S^a(t)`).
#' @noRd
compute_survival_if_matrix <- function(
  fit,
  intervention,
  times,
  prep,
  fit_idx,
  id_vec,
  unique_ids
) {
  pp_cf <- apply_intervention_pp(fit$pp_data, fit$treatment, intervention)

  ## Design matrix under the intervention. We strip the response side of the
  ## model's formula so model.matrix() accepts newdata even when the outcome
  ## column has been intervened upon in some weird downstream use case.
  tt <- stats::delete.response(stats::terms(fit$model))
  X_pp <- stats::model.matrix(tt, data = pp_cf)

  eta <- stats::predict(fit$model, newdata = pp_cf, type = "link")
  h <- stats::predict(fit$model, newdata = pp_cf, type = "response")
  mu_eta <- fit$model$family$mu.eta(eta)
  ## Per-row sensitivity of log(1 - h) wrt beta is -(s_row * X). For the
  ## logit link s_row simplifies to h itself; the general formula below
  ## handles any family providing mu_eta() and variance().
  s_row <- mu_eta / (1 - h)

  ## Within-id cumulative survival: at row (i, k), `.cf_surv = S_i(k)`.
  ## Within-id cumulative sensitivity: at row (i, k),
  ## `cum_SX[i, k, :] = sum_{l <= k} s_{i,l} * X_{i,l}`.
  id_col <- fit$id
  time_col <- fit$time
  data.table::setkeyv(pp_cf, c(id_col, time_col))
  ## Attach the per-row hazard as a column so the cumulative product inside
  ## `by = id_col` sees the grouped subset of `h` rather than the full
  ## enclosing vector (data.table's non-standard evaluation only scopes
  ## columns, not outer R variables, into the group expression).
  pp_cf[, .cf_hazard := h]
  pp_cf[, .cf_surv := cumprod(1 - .cf_hazard), by = c(id_col)]

  SX <- X_pp * s_row ## n_pp x p, row-scaled
  ## `ave()` per column gives a cumulative sum within id; applied column-wise
  ## this yields the n_pp x p cumulative sensitivity matrix.
  id_pp <- pp_cf[[id_col]]
  cum_SX <- apply(SX, 2L, function(v) stats::ave(v, id_pp, FUN = cumsum))

  ## For each user time t_j, pull the row per id at time == t_j and
  ## assemble the Ch1 and J_bar pieces.
  n_ids <- length(unique_ids)
  n_t <- length(times)
  p <- ncol(X_pp)

  Ch1_mat <- matrix(0, n_ids, n_t)
  J_bar_mat <- matrix(0, p, n_t)

  S_cf <- pp_cf[[".cf_surv"]]
  t_pp <- pp_cf[[time_col]]

  s_hat <- numeric(n_t)

  for (j in seq_len(n_t)) {
    rows_j <- which(t_pp == times[j])
    ## Guard against PP data that drops rows at certain times per id (Track
    ## A data is rectangular by construction but a defensive length check
    ## keeps a later chunk that allows ragged PP safe).
    if (length(rows_j) != n_ids) {
      rlang::abort(
        paste0(
          "survatr sandwich variance expected ",
          n_ids,
          " rows at time = ",
          times[j],
          ", got ",
          length(rows_j),
          "."
        ),
        class = "survatr_if_failed"
      )
    }
    S_at_tj <- S_cf[rows_j]
    s_hat[j] <- mean(S_at_tj)
    Ch1_mat[, j] <- S_at_tj - s_hat[j]

    ## A_j[i, :] = S_i(t_j) * cum_SX[rows_j, :][i, :]
    cum_SX_j <- cum_SX[rows_j, , drop = FALSE]
    A_j <- cum_SX_j * S_at_tj ## row-wise multiply
    J_bar_mat[, j] <- -colMeans(A_j)
  }

  ## Per-id beta IF: psi_per_id (n_ids x p) %*% B_inv.
  weighted_X <- prep$X_fit * prep$r_score ## n_fit x p
  id_fit <- id_vec[fit_idx]
  id_fit_f <- factor(id_fit, levels = unique_ids)
  ## `rowsum(..., reorder = FALSE)` keeps rows in the order of first
  ## occurrence in `id_fit_f`'s levels. We requested levels = unique_ids
  ## above so the result's row order matches `unique_ids` exactly. Ids
  ## with no fit rows (should be impossible under the chunk-1 risk-set
  ## rules but defend anyway) get a zero row.
  psi_per_id_raw <- rowsum(weighted_X, id_fit_f, reorder = FALSE)
  psi_per_id <- matrix(0, n_ids, p)
  rownames(psi_per_id) <- as.character(unique_ids)
  psi_per_id[rownames(psi_per_id_raw), ] <- psi_per_id_raw

  ## IF convention: theta_hat - theta ≈ (1/n) sum_i IF_i, so IF_i has
  ## magnitude O(1) and Var(theta_hat) ≈ (1/n^2) sum_i IF_i^2. For beta:
  ## beta_hat - beta ≈ -B_inv * sum_i psi_i = (1/n_ids) * sum_i (-n_ids *
  ## B_inv * psi_i), so the per-individual IF on beta is `-n_ids * B_inv
  ## * psi_i`. Without the n_ids scaling, crossprod(IF)/n_ids^2 would
  ## undersell the variance by a factor of n_ids^2.
  IF_beta_per_id <- n_ids * psi_per_id %*% prep$B_inv ## n_ids x p

  Ch2_mat <- IF_beta_per_id %*% J_bar_mat ## n_ids x n_t
  IF_mat <- Ch1_mat + Ch2_mat

  list(s_hat = s_hat, IF_mat = IF_mat)
}

#' Prepare shared pieces for the sandwich IF chain
#'
#' Rebuilds `fit_rows` / `fit_idx` from `fit$pp_data` (we stripped the
#' `.survatr_*` internal columns at the end of `surv_fit()`, so this is the
#' reverse), then calls `causatr:::prepare_model_if()` to get the bread
#' `B_inv` and working-residual vector `r_score`. Returned once and passed
#' in to `compute_survival_if_matrix()` for every intervention.
#'
#' @param fit A `survatr_fit`.
#'
#' @return A list with `prep` (the `causatr:::prepare_model_if()` output),
#'   `fit_idx`, `id_vec`, and `unique_ids`.
#' @noRd
prepare_sandwich_shared <- function(fit) {
  pp_work <- data.table::copy(fit$pp_data)
  fit_rows <- build_risk_set(
    data = pp_work,
    outcome = fit$outcome,
    id = fit$id,
    censoring = fit$censoring
  )
  fit_idx <- which(fit_rows)

  prep <- causatr:::prepare_model_if(
    model = fit$model,
    fit_idx = fit_idx,
    n_total = nrow(pp_work)
  )

  id_vec <- pp_work[[fit$id]]
  unique_ids <- unique(id_vec)

  list(
    prep = prep,
    fit_idx = fit_idx,
    id_vec = id_vec,
    unique_ids = unique_ids
  )
}
