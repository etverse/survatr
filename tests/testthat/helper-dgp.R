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

## Minimal three-period, five-id PP fixture used for row-level structural
## checks (risk-set construction, reserved-col guard, data-shape validation).
## id 1: event at t=2  -> at-risk rows: (1,1), (1,2). (1,3) dropped.
## id 2: no event.     -> at-risk rows: (2,1), (2,2), (2,3).
## id 3: censored t=2. -> at-risk rows: (3,1). (3,2), (3,3) dropped.
## id 4: event at t=1. -> at-risk rows: (4,1). (4,2), (4,3) dropped.
## id 5: no event.     -> at-risk rows: (5,1), (5,2), (5,3).
## Total at risk without censoring: 2+3+3+1+3 = 12.
## With censoring: 2+3+1+1+3 = 10.
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
