## Closed-form DGPs used across test files.

## Simulate a person-period dataset with constant discrete-time hazard
## `h` and no covariate or treatment effect. Used to pin `surv_fit()`'s
## family-switch oracle: the intercept of a pooled logistic fit with
## `time_formula = ~ 1` should converge to `qlogis(h)`.
##
## - n individuals, K periods (1..K).
## - Treatment A ~ Bernoulli(0.5), independent of outcome (no effect).
## - No covariates (confounders = ~ 1 at the call site).
## - Event indicator Y_k ~ Bernoulli(h) for rows still at risk.
## - Rows at/after the first event are present but marked Y_k = 0 and the
##   risk-set builder drops them.
## - No censoring.
sim_constant_hazard <- function(n = 2000L, K = 10L, h = 0.05, seed = 1L) {
  set.seed(seed)
  A_per_id <- stats::rbinom(n, 1L, 0.5)
  ## Draw a full K-period grid per id, then zero out Y after the first event.
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    Y <- stats::rbinom(K, 1L, h)
    first <- which(Y == 1L)[1L]
    if (!is.na(first) && first < K) {
      Y[(first + 1L):K] <- 0L
    }
    rows[[i]] <- data.table::data.table(
      id = i,
      t = seq_len(K),
      A = A_per_id[i],
      Y = Y
    )
  }
  data.table::rbindlist(rows)
}

## Confounded point-treatment survival DGP on a rectangular person-period grid.
## Used by the IPW tests: the treatment depends on a baseline confounder `L`
## that also drives the hazard, so a naive (unweighted) curve is biased and the
## IPW / gcomp adjusted curves should agree (and match the lmtp oracle).
##
## - n individuals, K periods (1..K).
## - Baseline confounder L ~ N(0, 1).
## - Treatment A ~ Bernoulli(plogis(gamma * L)) -- when `gamma = 0` the
##   treatment is independent of L (unconfounded), the stabilized weights
##   collapse to ~1, and the IPW curve matches the gcomp curve.
## - Discrete hazard logit h_{i,k} = alpha0 + alpha_t(k) + beta_A * A_i +
##   beta_L * L_i, with `alpha_t` small per-period offsets.
## - First-event truncation + rectangular padding (Y = 0 after the first event).
## - No censoring (keeps the lmtp oracle reshape simple).
sim_confounded_survival <- function(
  n = 2000L,
  K = 5L,
  seed = 1L,
  gamma = 0.8,
  alpha0 = -2.0,
  beta_A = -0.6,
  beta_L = 0.7
) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(gamma * L))
  ## Small deterministic per-period offsets so the baseline hazard is not flat.
  alpha_t <- seq(0, 0.3, length.out = K)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    eta <- alpha0 + alpha_t + beta_A * A[i] + beta_L * L[i]
    h <- stats::plogis(eta)
    Y <- stats::rbinom(K, 1L, h)
    first <- which(Y == 1L)[1L]
    if (!is.na(first) && first < K) {
      Y[(first + 1L):K] <- 0L
    }
    rows[[i]] <- data.table::data.table(
      id = i,
      t = seq_len(K),
      A = A[i],
      L = L[i],
      Y = Y
    )
  }
  data.table::rbindlist(rows)
}

## Analytical-ish truth for the confounded-survival DGP: the true marginal
## counterfactual survival S^a(t) = E_L[ prod_{k<=t} (1 - h_{k}(a, L)) ], obtained
## by Monte Carlo integration over L ~ N(0, 1) on a large fixed grid. The
## hazard form must match `sim_confounded_survival()`. Used to score sandwich CI
## coverage for the IPW estimator.
true_marginal_survival <- function(
  a,
  times,
  K = 5L,
  alpha0 = -2.0,
  beta_A = -0.6,
  beta_L = 0.7,
  n_int = 400000L,
  seed = 999L
) {
  set.seed(seed)
  L <- stats::rnorm(n_int)
  alpha_t <- seq(0, 0.3, length.out = K)
  ## h_mat[, k] = hazard at period k for each integration draw under A = a.
  h_mat <- vapply(
    seq_len(K),
    function(k) stats::plogis(alpha0 + alpha_t[k] + beta_A * a + beta_L * L),
    numeric(n_int)
  )
  surv_mat <- t(apply(1 - h_mat, 1L, cumprod)) # n_int x K cumulative product
  vapply(times, function(t) mean(surv_mat[, t]), numeric(1L))
}

## Minimal three-period, five-id PP fixture used for row-level structural
## checks (risk-set construction, reserved-col guard, data-shape validation).
## id 1: event at t=2  -> at-risk rows: (1,1), (1,2). (1,3) dropped.
## id 2: no event.     -> at-risk rows: (2,1), (2,2), (2,3).
## id 3: censored t=2. -> at-risk rows: (3,1). (3,2), (3,3) dropped.
## id 4: event at t=1. -> at-risk rows: (4,1). (4,2), (4,3) dropped.
## id 5: no event.     -> at-risk rows: (5,1), (5,2), (5,3).
## Total at risk without censoring: 2+3+3+1+3 = 12.
## With censoring: 2+3+1+1+3 = 10.

## DGP with INFORMATIVE censoring: censoring depends on L so the uncorrected
## row-filter IPW curve is biased. Setting `delta_cens = 0` gives
## non-informative censoring (IPCW curve == uncorrected curve).
##
## Model:
##   L ~ N(0, 1) (baseline confounder)
##   A ~ Bernoulli(logit^{-1}(gamma * L))    -- confounded treatment
##   h(t | A, L) = logit^{-1}(alpha0 + alpha_t * t + beta_A * A + beta_L * L)
##   P(C_t = 1 | A, L) = logit^{-1}(cens0 + delta_cens * L)  -- informative cens
##
## Returns a RECTANGULAR person-period data.table (padded with Y=0, C=0 after
## first event/censor) plus the true marginal survival for each arm
## (via Monte Carlo at n_int draws from L). True survival integrates over L,
## treating the censoring as non-informative in the DGP (the event process is
## defined unconditionally on censoring, and the IPCW estimator targets this
## marginal estimand).
sim_informative_censoring <- function(
  n = 2000L,
  K = 5L,
  seed = 1L,
  gamma = 0.7,
  alpha0 = -2.5,
  beta_A = -0.5,
  beta_L = 0.6,
  cens0 = -2.5,
  delta_cens = 0.8
) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(gamma * L))
  alpha_t <- seq(0, 0.4, length.out = K)

  rows <- vector("list", n)
  for (i in seq_len(n)) {
    eta_Y <- alpha0 + alpha_t + beta_A * A[i] + beta_L * L[i]
    h_Y <- stats::plogis(eta_Y)
    eta_C <- cens0 + delta_cens * L[i]
    h_C <- stats::plogis(eta_C)

    Y <- integer(K)
    C <- integer(K)
    alive <- TRUE
    for (k in seq_len(K)) {
      if (!alive) {
        break
      }
      Y[k] <- stats::rbinom(1L, 1L, h_Y[k])
      if (Y[k] == 1L) {
        alive <- FALSE
      } else {
        C[k] <- stats::rbinom(1L, 1L, h_C)
        if (C[k] == 1L) alive <- FALSE
      }
    }
    rows[[i]] <- data.table::data.table(
      id = i,
      t = seq_len(K),
      A = A[i],
      L = L[i],
      Y = Y,
      C = C
    )
  }
  data.table::rbindlist(rows)
}

## Marginal counterfactual survival for the informative-censoring DGP.
## Integrates the event process over L, ignoring censoring (target estimand).
true_marginal_survival_ipcw <- function(
  a,
  times,
  K = 5L,
  alpha0 = -2.5,
  beta_A = -0.5,
  beta_L = 0.6,
  n_int = 400000L,
  seed = 999L
) {
  set.seed(seed)
  L <- stats::rnorm(n_int)
  alpha_t <- seq(0, 0.4, length.out = K)
  h_mat <- vapply(
    seq_len(K),
    function(k) stats::plogis(alpha0 + alpha_t[k] + beta_A * a + beta_L * L),
    numeric(n_int)
  )
  surv_mat <- t(apply(1 - h_mat, 1L, cumprod))
  vapply(times, function(tt) mean(surv_mat[, tt]), numeric(1L))
}

fixture_small_pp <- function() {
  data.table::data.table(
    id = rep(1:5, each = 3L),
    t = rep(1:3, times = 5L),
    A = rep(c(1L, 0L, 1L, 0L, 1L), each = 3L),
    L = rep(c(0.1, -0.3, 0.2, 0.4, -0.1), each = 3L),
    Y = c(
      0L,
      1L,
      0L, ## id 1
      0L,
      0L,
      0L, ## id 2
      0L,
      0L,
      0L, ## id 3 (censored at t=2 per cens column)
      1L,
      0L,
      0L, ## id 4
      0L,
      0L,
      0L ## id 5
    ),
    cens = c(
      0L,
      0L,
      0L, ## id 1
      0L,
      0L,
      0L, ## id 2
      0L,
      1L,
      0L, ## id 3 censored at t=2
      0L,
      0L,
      0L, ## id 4
      0L,
      0L,
      0L ## id 5
    )
  )
}
