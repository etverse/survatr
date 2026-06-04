## Oracles + DGP for Track B (ICE-hazard longitudinal survival) tests.
##
## The DGP is a treatment-confounder-feedback discrete-time survival process
## (Hernan & Robins Ch. 21 / Table 20.1 structure): a time-varying confounder
## L_k drives both the time-varying treatment A_k and the hazard, and past
## treatment A_{k-1} drives future L_k -- so a naive (unadjusted) curve is
## biased and only the longitudinal g-formula (ICE) recovers the truth.
##
## The PRIMARY oracle is the forward-simulation Monte-Carlo g-formula truth
## (`true_ice_survival()`): under a static strategy A_k == a* it forward-
## simulates the intervened world (treatment a* at every period, with L_k still
## evolving under a*) and averages the cumulative survival. This is exactly the
## g-formula functional the ICE estimator targets, so it is an analytical
## reference, not a smoke test. `oracle_lmtp_survival()` is a secondary
## external cross-check (skipped when lmtp is absent).

## Default DGP parameters (shared by the simulator and the truth so they stay
## consistent). `a_prev0` is the pre-study treatment A_0 (the lag at period 1).
.ice_dgp_params <- function() {
  list(
    g_L = 0.6, # effect of current L on current treatment
    g_A = -0.3, # effect of lagged treatment on current treatment
    mu_L = 0.4, # autoregression of the confounder
    mu_A = 0.5, # feedback of lagged treatment onto the confounder
    a0 = -2.0, # hazard intercept
    bL = 0.8, # effect of current L on the hazard
    bA = -0.7, # protective effect of current treatment on the hazard
    a_prev0 = 0
  )
}

## Simulate an OBSERVATIONAL person-period dataset with treatment-confounder
## feedback. Rectangular grid (1..K) per id; rows at/after the first event are
## padded (Y = 0, covariates carried forward) so the risk-set builder drops
## them from the fit while keeping the grid rectangular.
sim_ice_feedback <- function(
  n = 800L,
  K = 4L,
  seed = 1L,
  p = .ice_dgp_params()
) {
  set.seed(seed)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    L <- numeric(K)
    A <- numeric(K)
    Y <- integer(K)
    l_prev <- 0
    a_prev <- p$a_prev0
    failed_at <- NA_integer_
    for (k in seq_len(K)) {
      L[k] <- stats::rnorm(1, p$mu_L * l_prev + p$mu_A * a_prev, 1)
      A[k] <- stats::rbinom(1, 1L, stats::plogis(p$g_L * L[k] + p$g_A * a_prev))
      h <- stats::plogis(p$a0 + p$bL * L[k] + p$bA * A[k])
      Y[k] <- stats::rbinom(1, 1L, h)
      l_prev <- L[k]
      a_prev <- A[k]
      if (Y[k] == 1L) {
        failed_at <- k
        break
      }
    }
    ## Pad post-event rows (rectangular grid). Carry covariates / treatment
    ## forward; these rows are never at risk so their values don't enter any
    ## fit or prediction, but must be non-NA (predictor NA is rejected).
    if (!is.na(failed_at) && failed_at < K) {
      for (j in (failed_at + 1L):K) {
        L[j] <- L[failed_at]
        A[j] <- A[failed_at]
        Y[j] <- 0L
      }
    }
    rows[[i]] <- data.table::data.table(
      id = i,
      t = seq_len(K),
      L = L,
      A = A,
      Y = Y
    )
  }
  data.table::rbindlist(rows)
}

## Forward-simulation Monte-Carlo g-formula truth for a STATIC strategy
## A_k == a at every period. Simulates the intervened world (treatment fixed
## at `a`, L_k evolving under `a`) for a large sample and returns the marginal
## counterfactual survival S^a(t) = E[ prod_{k<=t} (1 - h_k) ] at each `times`.
true_ice_survival <- function(
  a,
  times,
  K = 4L,
  p = .ice_dgp_params(),
  n_int = 200000L,
  seed = 7L
) {
  set.seed(seed)
  surv <- matrix(1, n_int, K)
  l_prev <- rep(0, n_int)
  a_prev <- rep(p$a_prev0, n_int) # pre-study A_0
  s_cum <- rep(1, n_int)
  for (k in seq_len(K)) {
    l_k <- stats::rnorm(n_int, p$mu_L * l_prev + p$mu_A * a_prev, 1)
    h_k <- stats::plogis(p$a0 + p$bL * l_k + p$bA * a)
    s_cum <- s_cum * (1 - h_k)
    surv[, k] <- s_cum
    l_prev <- l_k
    a_prev <- a # treatment forced to `a` at every period
  }
  vapply(times, function(t) mean(surv[, t]), numeric(1L))
}

## Forward-simulation g-formula truth for a DYNAMIC strategy "treat iff the
## current confounder exceeds `thresh`" (A_k = 1{L_k > thresh}). A genuinely
## longitudinal, non-static strategy: treatment responds to the evolving
## covariate, which in turn responds to past treatment (feedback). Simulates the
## intervened world and returns S^d(t) at each `times`.
true_ice_dynamic <- function(
  thresh,
  times,
  K = 4L,
  p = .ice_dgp_params(),
  n_int = 200000L,
  seed = 8L
) {
  set.seed(seed)
  surv <- matrix(1, n_int, K)
  l_prev <- rep(0, n_int)
  a_prev <- rep(p$a_prev0, n_int)
  s_cum <- rep(1, n_int)
  for (k in seq_len(K)) {
    l_k <- stats::rnorm(n_int, p$mu_L * l_prev + p$mu_A * a_prev, 1)
    a_k <- as.numeric(l_k > thresh) # the dynamic rule, applied to current L
    h_k <- stats::plogis(p$a0 + p$bL * l_k + p$bA * a_k)
    s_cum <- s_cum * (1 - h_k)
    surv[, k] <- s_cum
    l_prev <- l_k
    a_prev <- a_k
  }
  vapply(times, function(t) mean(surv[, t]), numeric(1L))
}

## Secondary external oracle: lmtp::lmtp_tmle(outcome_type = "survival") for a
## STATIC longitudinal strategy A_k == a, with a TIME-VARYING confounder.
## Mirrors the existing IPW lmtp oracle (helper-ipw-survival-oracle.R): wide
## reshape, cumulative-event outcome nodes, always-observed censoring nodes,
## one fit per horizon, and the S7 `ife` estimate at `@x`. Adapted for the
## longitudinal case: per-period treatment nodes A_1..A_tau and per-period
## time-varying covariate L_k. Returns the S^a(t) vector over `times` (default
## 2:K; lmtp's survival fit needs >= 2 outcome nodes, so t = 1 is left to the
## forward-sim truth), or NULL if lmtp is unavailable / any horizon fails.
## Validated against `true_ice_survival()` at test design time before being
## pinned as an oracle.
oracle_lmtp_survival <- function(dat, a, K = 4L, times = NULL) {
  if (!requireNamespace("lmtp", quietly = TRUE)) {
    return(NULL)
  }
  if (is.null(times)) {
    times <- seq.int(2L, K)
  }
  dt <- data.table::as.data.table(dat)
  ## Wide: one row per id with A_k, L_k, Y_k for each period.
  w <- data.table::dcast(dt, id ~ t, value.var = c("A", "L", "Y"))
  ## Survival outcome: cumulative (monotone) event indicator.
  for (k in 2:K) {
    w[[paste0("Y_", k)]] <- pmax(
      w[[paste0("Y_", k - 1L)]],
      w[[paste0("Y_", k)]]
    )
  }
  ## Always-observed censoring nodes (the DGP has no censoring).
  for (k in seq_len(K)) {
    w[[paste0("C_", k)]] <- 1L
  }
  wide_df <- as.data.frame(w)

  out <- vapply(
    times,
    function(tau) {
      fit_lmtp <- tryCatch(
        lmtp::lmtp_tmle(
          data = wide_df,
          trt = paste0("A_", seq_len(tau)),
          outcome = paste0("Y_", seq_len(tau)),
          time_vary = lapply(seq_len(tau), function(k) paste0("L_", k)),
          cens = paste0("C_", seq_len(tau)),
          shift = function(data, trt) rep(a, nrow(data)),
          outcome_type = "survival",
          folds = 1L
        ),
        error = function(e) NULL
      )
      if (is.null(fit_lmtp)) {
        return(NA_real_)
      }
      ## S7 `ife` object: survival point estimate at horizon `tau` is `@x`.
      as.numeric(fit_lmtp$estimate@x)
    },
    numeric(1L)
  )
  if (anyNA(out)) {
    return(NULL)
  }
  out
}
