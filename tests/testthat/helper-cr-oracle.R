## Competing-risks oracles and DGP shared across the CR test files.

## Two-cause discrete-time competing-risks DGP on a rectangular person-period
## grid. At each at-risk period the per-period outcome is a 3-category
## multinomial: no event with probability 1 - h1 - h2, a cause-1 event with
## probability h1, a cause-2 event with probability h2. A logistic regression of
## 1{event == j} on the at-risk rows therefore estimates the cause-specific
## hazard h_j exactly, so a pooled-logistic cause-specific model is correctly
## specified and the closed-form CIF below is the truth (up to MC error).
##
## - n individuals, K periods (1..K).
## - Baseline confounder L ~ N(0, 1); treatment A ~ Bernoulli(plogis(gamma*L)).
## - Cause-j hazard logit: qlogis(h_j) + beta_jA * A + beta_jL * L.
## - With beta_*A = beta_*L = gamma = 0 the hazards are constant and treatment-
##   and covariate-free: the analytic closed form `analytic_cif()` applies and
##   both arms are identical (no-effect DGP for coverage of 0 / 1).
## - First-event truncation + rectangular padding (event = 0 after the first
##   event of any cause). No administrative censoring.
sim_two_cause_constant_hazard <- function(
  n = 5000L,
  K = 8L,
  h1 = 0.06,
  h2 = 0.04,
  seed = 1L,
  beta1_A = 0,
  beta2_A = 0,
  beta1_L = 0,
  beta2_L = 0,
  gamma = 0
) {
  set.seed(seed)
  L <- stats::rnorm(n)
  A <- stats::rbinom(n, 1L, stats::plogis(gamma * L))
  l1 <- stats::qlogis(h1)
  l2 <- stats::qlogis(h2)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    p1 <- stats::plogis(l1 + beta1_A * A[i] + beta1_L * L[i])
    p2 <- stats::plogis(l2 + beta2_A * A[i] + beta2_L * L[i])
    event <- integer(K)
    done <- FALSE
    for (k in seq_len(K)) {
      if (done) {
        event[k] <- 0L
        next
      }
      ## Multinomial draw: (no event, cause 1, cause 2). When the conditional
      ## hazards are large enough that p1 + p2 > 1 we cap "no event" at 0; the
      ## constant-hazard oracle uses small hazards so this never binds there.
      u <- stats::runif(1L)
      if (u < p1) {
        event[k] <- 1L
        done <- TRUE
      } else if (u < p1 + p2) {
        event[k] <- 2L
        done <- TRUE
      } else {
        event[k] <- 0L
      }
    }
    rows[[i]] <- data.table::data.table(
      id = i,
      t = seq_len(K),
      A = A[i],
      L = L[i],
      event = event
    )
  }
  data.table::rbindlist(rows)
}

## Closed-form CIF / survival for the constant-hazard two-cause DGP. With
## summed hazard H = h1 + h2 and S(t) = (1 - H)^t,
##   F_j(t) = sum_{k<=t} S(k-1) h_j = (h_j / H) * (1 - (1 - H)^t).
## The identity F1(t) + F2(t) + S(t) = 1 is exact.
analytic_cr <- function(times, h1, h2) {
  hsum <- h1 + h2
  surv <- (1 - hsum)^times
  data.frame(
    time = times,
    F1 = (h1 / hsum) * (1 - surv),
    F2 = (h2 / hsum) * (1 - surv),
    S = surv
  )
}

## riskRegression::ate g-formula CIF oracle. Collapses person-period data to
## one-row-per-id, fits cause-specific Cox models (CSC) for all causes, and
## returns the standardised CIF under each treatment arm via ate(). Requires
## the riskRegression package (Suggests). Treatment is coerced to factor
## because ate() enforces this.
##
## Returns a data.table with columns: cause, time, treatment (factor), cif.
## `treatment` levels match the original 0/1 coding cast as character ("0"/"1").
rr_ate_cif_oracle <- function(
  pp,
  times,
  id = "id",
  time = "t",
  event = "event",
  treatment = "A",
  confounders = "L"
) {
  dt <- data.table::as.data.table(pp)
  data.table::setkeyv(dt, c(id, time))
  K <- max(dt[[time]])
  ## Collapse PP to one row per id: first-event time + cause (0 = censored).
  per_id <- dt[,
    {
      ev <- which(.SD[[1L]] != 0L)
      if (length(ev) > 0L) {
        list(stop = .SD[[2L]][ev[1L]], event_type = .SD[[1L]][ev[1L]])
      } else {
        list(stop = K, event_type = 0L)
      }
    },
    by = c(id),
    .SDcols = c(event, time)
  ]
  ## Merge back baseline covariates (first row per id).
  baseline <- dt[, .SD[1L], by = id, .SDcols = c(treatment, confounders)]
  per_id <- merge(per_id, baseline, by = id)
  ## ate() requires treatment to be a factor.
  per_id[[treatment]] <- factor(per_id[[treatment]])
  rhs <- paste(c(treatment, confounders), collapse = " + ")
  ## Build formula with Hist resolved via an explicit env so the bare name
  ## "Hist" in the call satisfies riskRegression:::SurvResponseVar's pattern
  ## check while avoiding a namespace-qualified LHS that breaks the check.
  fml_env <- new.env(parent = baseenv())
  fml_env$Hist <- prodlim::Hist
  fml <- stats::as.formula(
    paste0("Hist(stop, event_type) ~ ", rhs),
    env = fml_env
  )
  csc_fit <- riskRegression::CSC(fml, data = per_id)
  ## ate() re-evaluates csc_fit$call$formula; store the object so the lookup
  ## succeeds outside the function scope.
  csc_fit$call$formula <- fml
  causes <- sort(unique(per_id$event_type[per_id$event_type != 0L]))
  out <- lapply(causes, function(j) {
    res <- riskRegression::ate(
      csc_fit,
      data = per_id,
      treatment = treatment,
      times = times,
      cause = j,
      se = FALSE,
      verbose = FALSE
    )
    mr <- data.table::as.data.table(res$meanRisk)
    mr[, cause := j]
    mr[, .(cause, time, treatment, cif = estimate)]
  })
  data.table::rbindlist(out)
}

## Aalen-Johansen CIF oracle via `survival::survfit` on the collapsed
## (one-row-per-id) data. Returns a long data.frame `cause | time | cif` for the
## requested integer event times. Used as an empirical (non-closed-form)
## cross-check against the saturated (`~ factor(t)`) gcomp CIF, which reproduces
## the per-period empirical cause-specific hazards AJ uses. Requires `survival`.
aj_cif_oracle <- function(pp, times, id = "id", time = "t", event = "event") {
  ## Collapse PP => one row per id: event period (or K, censored) + status.
  dt <- data.table::as.data.table(pp)
  data.table::setkeyv(dt, c(id, time))
  K <- max(dt[[time]])
  per_id <- dt[,
    {
      ev <- which(.SD[[1L]] != 0)
      if (length(ev) > 0L) {
        list(stop = .SD[[2L]][ev[1L]], status = .SD[[1L]][ev[1L]])
      } else {
        list(stop = K, status = 0L)
      }
    },
    by = c(id),
    .SDcols = c(event, time)
  ]
  fit <- survival::survfit(
    survival::Surv(stop, factor(status)) ~ 1,
    data = per_id
  )
  ## `pstate` columns are the states in `fit$states`: the event-free state
  ## (labelled "(s0)" / "" by survfitms) plus one per cause ("1", "2", ...). The
  ## cause-j CIF is the probability in state "j"; keep only the integer-labelled
  ## (event) states.
  sm <- summary(fit, times = times, extend = TRUE)
  states <- fit$states
  pst <- sm$pstate
  colnames(pst) <- states
  event_states <- states[grepl("^[0-9]+$", states) & states != "0"]
  out <- lapply(event_states, function(s) {
    data.frame(
      cause = as.integer(s),
      time = sm$time,
      cif = pst[, s]
    )
  })
  do.call(rbind, out)
}
