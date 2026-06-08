#' Shared sandwich pieces for the competing-risks IF chain
#'
#' Rebuilds the shared all-cause risk set from `fit$pp_data` (the `.survatr_*`
#' columns were stripped at the end of `surv_fit()`), then calls
#' `causatr:::prepare_model_if()` **once per cause** to get the per-cause bread
#' `B_inv` and working-residual score. The bread is block-diagonal across
#' causes (each cause model's score involves only its own coefficients), so each
#' cause's `prepare_model_if()` is independent; the cross-cause correlation
#' enters later when the per-individual IFs (summed over causes) cross-multiply.
#'
#' @param fit A competing-risks `survatr_fit`.
#'
#' @returns A list with `prep_by_cause` (named list of `prepare_model_if()`
#'   outputs, one per cause), `fit_idx` (`which(fit_rows)`), `id_vec` (id per PP
#'   row), and `unique_ids`.
#' @family competing-risks
#' @noRd
prepare_cr_sandwich_shared <- function(fit) {
  pp_work <- data.table::copy(fit$pp_data)
  ## All-cause event indicator => shared risk set (same as the fit step).
  pp_work[, .survatr_any_event := as.integer(pp_work[[fit$competing]] != 0)]
  fit_rows <- build_risk_set(
    data = pp_work,
    outcome = ".survatr_any_event",
    id = fit$id,
    censoring = fit$censoring
  )
  fit_idx <- which(fit_rows)

  prep_by_cause <- lapply(fit$causes, function(j) {
    m <- fit$cause_models[[as.character(j)]]
    ## NA rows are rejected at fit time; assert the model dropped none so the
    ## fit-row alignment used below holds (mirrors the single-event guard).
    if (length(m$na.action) > 0L) {
      rlang::abort(
        "Competing-risks sandwich: a cause model has non-empty na.action.",
        class = "survatr_if_failed"
      )
    }
    causatr:::prepare_model_if(
      model = m,
      fit_idx = fit_idx,
      n_total = nrow(pp_work)
    )
  })
  names(prep_by_cause) <- as.character(fit$causes)

  list(
    prep_by_cause = prep_by_cause,
    fit_idx = fit_idx,
    id_vec = pp_work[[fit$id]],
    unique_ids = unique(pp_work[[fit$id]])
  )
}

#' Per-intervention influence-function building blocks (competing risks)
#'
#' Precomputes everything the per-cause CIF IF and the all-cause survival IF
#' share for one intervention: the counterfactual design matrix, per-cause
#' hazards / `mu.eta`, the all-cause survival and its lag, the per-cause
#' cumulative sensitivities, and the per-individual coefficient IFs. Computing
#' these once (rather than per target cause) keeps the `J`-cause sandwich
#' tractable.
#'
#' The per-cause sensitivity is `s^(j')_l = mu_eta^(j')_l / (1 - H_l)` with the
#' **all-cause** denominator `1 - H` (`H = sum_j h^(j)`); unlike the
#' single-event case (chunk 3) the `(1 - h)` does not cancel `mu_eta = h(1 - h)`
#' because the survival uses the summed hazard. `cumSX^(j')` is the within-id
#' cumulative `sum_{l<=k} s^(j')_l X_l`; `cumSX_lag^(j') = cumSX^(j') - SX^(j')`
#' is the same summed only up to `k-1` (needed by the CIF derivative, which
#' weights by `S(k-1)`).
#'
#' @param fit A competing-risks `survatr_fit`.
#' @param intervention A single `causatr_intervention`.
#' @param shared Output of `prepare_cr_sandwich_shared()`.
#'
#' @returns A list of building blocks: `id_pp`, `t_pp`, `surv`, `surv_lag`,
#'   `h` (named list of per-cause hazard vectors), `mu` (per-cause `mu.eta`),
#'   `cumSX` / `cumSX_lag` (named lists of `n_pp x p` matrices), `IF_beta`
#'   (named list of `n_ids x p` per-individual coefficient IFs), `X_pp`,
#'   `causes`, `n_ids`.
#' @family competing-risks
#' @noRd
cr_intervention_if_pieces <- function(fit, intervention, shared) {
  id_col <- fit$id
  time_col <- fit$time
  causes <- fit$causes

  pp_cf <- apply_intervention_pp(fit$pp_data, fit$treatment, intervention)
  pp_cf <- cr_augment_pp(pp_cf, fit$cause_models, causes, id_col, time_col)

  id_pp <- pp_cf[[id_col]]
  t_pp <- pp_cf[[time_col]]
  surv <- pp_cf[[".cf_surv"]]
  surv_lag <- pp_cf[[".cf_surv_lag"]]
  all_haz <- pp_cf[[".cf_all_haz"]]

  ## All cause models share one design basis; build it once (lpmatrix for a
  ## gam, terms design for a glm) so it is conformable with every cause's bread.
  X_pp <- causatr:::iv_design_matrix(fit$cause_models[[1L]], pp_cf)

  h <- list()
  mu <- list()
  cumSX <- list()
  cumSX_lag <- list()
  for (j in causes) {
    jc <- as.character(j)
    m <- fit$cause_models[[jc]]
    eta <- as.numeric(stats::predict(m, newdata = pp_cf, type = "link"))
    h[[jc]] <- pp_cf[[paste0(".cf_haz_", j)]]
    mu[[jc]] <- m$family$mu.eta(eta)
    ## s^(j')_row = mu_eta^(j') / (1 - H) -- all-cause denominator.
    s_row <- mu[[jc]] / (1 - all_haz)
    SX <- X_pp * s_row
    cum <- apply(SX, 2L, function(v) stats::ave(v, id_pp, FUN = cumsum))
    cumSX[[jc]] <- cum
    cumSX_lag[[jc]] <- cum - SX
  }

  ## Per-individual coefficient IFs, one block per cause (block-diagonal bread).
  unique_ids <- shared$unique_ids
  n_ids <- length(unique_ids)
  p <- ncol(X_pp)
  id_fit <- shared$id_vec[shared$fit_idx]
  id_fit_f <- factor(id_fit, levels = unique_ids)
  IF_beta <- list()
  for (j in causes) {
    jc <- as.character(j)
    prep <- shared$prep_by_cause[[jc]]
    weighted_X <- prep$X_fit * prep$r_score
    psi_raw <- rowsum(weighted_X, id_fit_f, reorder = FALSE)
    psi <- matrix(0, n_ids, p)
    rownames(psi) <- as.character(unique_ids)
    psi[rownames(psi_raw), ] <- psi_raw
    ## IF on beta: beta_hat - beta ~ (1/n_ids) sum_i (-n_ids B_inv psi_i), so
    ## the per-individual IF carries the n_ids factor (see chunk 3).
    IF_beta[[jc]] <- n_ids * psi %*% prep$B_inv
  }

  list(
    id_pp = id_pp,
    t_pp = t_pp,
    surv = surv,
    surv_lag = surv_lag,
    h = h,
    mu = mu,
    cumSX = cumSX,
    cumSX_lag = cumSX_lag,
    IF_beta = IF_beta,
    X_pp = X_pp,
    causes = causes,
    n_ids = n_ids
  )
}

#' Per-individual IF matrix for the cause-`j` CIF `F^(j),a(t)`
#'
#' Implements the stacked cross-time delta chain for the cumulative incidence
#' `F^(j)_i(t) = sum_{k<=t} S_i(k-1) h^(j)_{i,k}`:
#'
#' ```
#' IF_i(t) = (F^(j)_i(t) - F_bar(t))                          ## Ch1
#'         + sum_{j'} IF_beta^(j')_i %*% J_bar^(j')(t)        ## Ch2
#' ```
#'
#' with the per-cause population derivative `J_bar^(j')(t) = mean_i d F^(j)_i(t) /
#' d beta_{j'}` assembled from
#' `d F^(j)_i(t)/d beta_{j'} = sum_{k<=t} [ -S(k-1) cumSX^(j')(k-1) h^(j)_k +
#' 1{j'=j} S(k-1) mu^(j)_k X_k ]`.
#'
#' @param pieces Output of `cr_intervention_if_pieces()`.
#' @param cause_j Integer target cause label.
#' @param times User time grid (sorted unique).
#'
#' @returns A list with `f_hat` (length `|times|`) and `IF_mat`
#'   (`n_ids x |times|`).
#' @family competing-risks
#' @noRd
compute_cif_if_matrix <- function(pieces, cause_j, times) {
  jc <- as.character(cause_j)
  id_pp <- pieces$id_pp
  t_pp <- pieces$t_pp
  surv_lag <- pieces$surv_lag
  h_j <- pieces$h[[jc]]
  mu_j <- pieces$mu[[jc]]
  X_pp <- pieces$X_pp
  n_ids <- pieces$n_ids
  n_t <- length(times)
  p <- ncol(X_pp)

  ## Per-individual CIF F^(j)_i(k): cumulative within id of S(k-1) h^(j)_k.
  f_inc <- surv_lag * h_j
  cum_f <- stats::ave(f_inc, id_pp, FUN = cumsum)

  ## Cause-j derivative weights shared across all j': the term_A scalar
  ## -S(k-1) h^(j)_k (row-multiplies cumSX_lag^(j')) and the term_B scalar
  ## S(k-1) mu^(j)_k (row-multiplies X, only when j' == j).
  wA <- -(surv_lag * h_j)
  wB <- surv_lag * mu_j

  ## J_bar^(j')(t) for every derivative cause j', accumulated into Ch2.
  Ch2 <- matrix(0, n_ids, n_t)
  for (jp in pieces$causes) {
    jpc <- as.character(jp)
    ## D^(j')_{i,k} = wA * cumSX_lag^(j') + 1{j'=j} wB * X  (n_pp x p).
    D <- pieces$cumSX_lag[[jpc]] * wA
    if (identical(jp, cause_j)) {
      D <- D + X_pp * wB
    }
    ## Cumulate within id over k, then pull rows at each requested time.
    cum_D <- apply(D, 2L, function(v) stats::ave(v, id_pp, FUN = cumsum))
    J_bar <- matrix(0, p, n_t)
    for (tj in seq_len(n_t)) {
      rows_tj <- which(t_pp == times[tj])
      J_bar[, tj] <- colMeans(cum_D[rows_tj, , drop = FALSE])
    }
    Ch2 <- Ch2 + pieces$IF_beta[[jpc]] %*% J_bar
  }

  ## Ch1 and the point estimate.
  Ch1 <- matrix(0, n_ids, n_t)
  f_hat <- numeric(n_t)
  for (tj in seq_len(n_t)) {
    rows_tj <- which(t_pp == times[tj])
    f_at <- cum_f[rows_tj]
    f_hat[tj] <- mean(f_at)
    Ch1[, tj] <- f_at - f_hat[tj]
  }

  list(f_hat = f_hat, IF_mat = Ch1 + Ch2)
}

#' Per-individual IF matrix for all-cause survival `S^a(t)` (competing risks)
#'
#' All-cause survival from the summed cause-specific hazards. The chain is the
#' chunk-3 single-event form generalised to multiple causes:
#' `d S_i(t)/d beta_{j'} = -S_i(t) cumSX^(j')_i(t)`, so
#' `J_bar_S^(j')(t) = -mean_i [ S_i(t) cumSX^(j')_i(t) ]`. Reduces to the
#' single-event survival IF when `J = 1`.
#'
#' @param pieces Output of `cr_intervention_if_pieces()`.
#' @param times User time grid.
#'
#' @returns A list with `s_hat` (length `|times|`) and `IF_mat`
#'   (`n_ids x |times|`).
#' @family competing-risks
#' @noRd
compute_cr_survival_if_matrix <- function(pieces, times) {
  id_pp <- pieces$id_pp
  t_pp <- pieces$t_pp
  surv <- pieces$surv
  n_ids <- pieces$n_ids
  n_t <- length(times)
  p <- ncol(pieces$X_pp)

  Ch1 <- matrix(0, n_ids, n_t)
  s_hat <- numeric(n_t)
  ## J_bar_S^(j')(t) per derivative cause, accumulated.
  J_bar <- lapply(pieces$causes, function(.) matrix(0, p, n_t))
  names(J_bar) <- as.character(pieces$causes)

  for (tj in seq_len(n_t)) {
    rows_tj <- which(t_pp == times[tj])
    s_at <- surv[rows_tj]
    s_hat[tj] <- mean(s_at)
    Ch1[, tj] <- s_at - s_hat[tj]
    for (jp in pieces$causes) {
      jpc <- as.character(jp)
      cum_at <- pieces$cumSX[[jpc]][rows_tj, , drop = FALSE]
      ## A_j[i, ] = S_i(t) * cumSX^(j')_i(t); J_bar = -colMeans(A_j).
      J_bar[[jpc]][, tj] <- -colMeans(cum_at * s_at)
    }
  }

  Ch2 <- matrix(0, n_ids, n_t)
  for (jp in pieces$causes) {
    jpc <- as.character(jp)
    Ch2 <- Ch2 + pieces$IF_beta[[jpc]] %*% J_bar[[jpc]]
  }

  list(s_hat = s_hat, IF_mat = Ch1 + Ch2)
}

#' Fill sandwich SEs / CIs into a competing-risks `survatr_result`
#'
#' Builds the per-(intervention, cause) CIF influence-function matrices (or the
#' per-intervention all-cause survival IF), aggregates them into pointwise SEs
#' and Wald CIs at `conf_level`, and replaces the `NA_real_` placeholders in
#' `estimates` / `contrasts`. The CIF-difference contrast IF is
#' `IF_F_a1 - IF_F_ref`; the CIF-ratio is built on the log scale and
#' exponentiated (strictly positive CI), mirroring the single-event
#' `risk_ratio`.
#'
#' @param fit A competing-risks `survatr_fit`.
#' @param estimates,contrasts The point-estimate tables (with a `cause` column).
#' @param interventions Named list of interventions.
#' @param times User time grid.
#' @param type Competing-risks contrast type.
#' @param reference Reference intervention name, or `NULL`.
#' @param causes Integer vector of reported causes.
#' @param conf_level Scalar in (0, 1).
#' @param shared Output of `prepare_cr_sandwich_shared()`.
#'
#' @returns A list `list(estimates, contrasts)` with SE / CI columns filled.
#' @family competing-risks
#' @noRd
fill_sandwich_ses_cr <- function(
  fit,
  estimates,
  contrasts,
  interventions,
  times,
  type,
  reference,
  causes,
  conf_level,
  shared
) {
  estimates <- data.table::copy(estimates)
  contrasts <- data.table::copy(contrasts)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  n_ids <- length(shared$unique_ids)
  iv_names <- names(interventions)

  ## Precompute per-intervention building blocks once.
  pieces <- lapply(iv_names, function(iv) {
    cr_intervention_if_pieces(fit, interventions[[iv]], shared)
  })
  names(pieces) <- iv_names

  if (type %in% c("survival", "risk")) {
    ## All-cause estimands: one IF per intervention; no pairwise contrasts.
    for (iv in iv_names) {
      ifm <- compute_cr_survival_if_matrix(pieces[[iv]], times)
      se_vec <- sqrt(pmax(diag(crossprod(ifm$IF_mat)) / n_ids^2, 0))
      tgt <- if (type == "survival") "s_hat" else "risk_hat"
      pt <- estimates[get("intervention") == iv & is.na(get("cause")), get(tgt)]
      estimates[
        get("intervention") == iv & is.na(get("cause")),
        `:=`(
          se = se_vec,
          ci_lower = pt - z * se_vec,
          ci_upper = pt + z * se_vec
        )
      ]
    }
    return(list(estimates = estimates, contrasts = contrasts))
  }

  ## CIF estimands: one IF per (intervention, cause).
  if_by <- lapply(iv_names, function(iv) {
    per_cause <- lapply(causes, function(j) {
      compute_cif_if_matrix(pieces[[iv]], j, times)
    })
    names(per_cause) <- as.character(causes)
    per_cause
  })
  names(if_by) <- iv_names

  ## Per-intervention per-cause SEs on the estimates table.
  for (iv in iv_names) {
    for (j in causes) {
      ifm <- if_by[[iv]][[as.character(j)]]
      se_vec <- sqrt(pmax(diag(crossprod(ifm$IF_mat)) / n_ids^2, 0))
      pt <- estimates[
        get("intervention") == iv & get("cause") == j,
        get("cif_hat")
      ]
      estimates[
        get("intervention") == iv & get("cause") == j,
        `:=`(
          se = se_vec,
          ci_lower = pt - z * se_vec,
          ci_upper = pt + z * se_vec
        )
      ]
    }
  }

  if (type == "cif" || nrow(contrasts) == 0L) {
    return(list(estimates = estimates, contrasts = contrasts))
  }

  ## Pairwise CIF contrasts.
  other_names <- setdiff(iv_names, reference)
  for (a1 in other_names) {
    cn <- paste0(a1, " vs ", reference)
    for (j in causes) {
      jc <- as.character(j)
      if_a1 <- if_by[[a1]][[jc]]$IF_mat
      if_ref <- if_by[[reference]][[jc]]$IF_mat
      f_a1 <- if_by[[a1]][[jc]]$f_hat
      f_ref <- if_by[[reference]][[jc]]$f_hat
      sel <- contrasts[, get("contrast") == cn & get("cause") == j]
      if (type == "cif_difference") {
        if_diff <- if_a1 - if_ref
        se_vec <- sqrt(pmax(diag(crossprod(if_diff)) / n_ids^2, 0))
        est_vec <- contrasts[sel, get("estimate")]
        contrasts[
          sel,
          `:=`(
            se = se_vec,
            ci_lower = est_vec - z * se_vec,
            ci_upper = est_vec + z * se_vec
          )
        ]
      } else {
        ## cif_ratio: delta on the log scale.
        ##   IF_{log RR} = (1/F_a1) IF_a1 - (1/F_ref) IF_ref
        if_log <- sweep(if_a1, 2L, f_a1, "/") -
          sweep(if_ref, 2L, f_ref, "/")
        if (any(f_a1 == 0) || any(f_ref == 0)) {
          bad <- which(f_a1 == 0 | f_ref == 0)
          if_log[, bad] <- NA_real_
        }
        se_log <- sqrt(pmax(diag(crossprod(if_log)) / n_ids^2, 0))
        ratio_vec <- contrasts[sel, get("estimate")]
        contrasts[
          sel,
          `:=`(
            se = ratio_vec * se_log,
            ci_lower = exp(log(ratio_vec) - z * se_log),
            ci_upper = exp(log(ratio_vec) + z * se_log)
          )
        ]
      }
    }
  }

  list(estimates = estimates, contrasts = contrasts)
}
