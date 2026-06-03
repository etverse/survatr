# Chunk 18 — Recurrent events + multi-state models

> **Status: ⬜ Not started · LARGE LIFT (data-model change)**
> **Depends on:** Chunk 2 (contrast spine), Chunk 3 (sandwich IF), Chunk 7
> (cause-specific / CIF machinery generalizes to transition-specific hazards).
> **Oracle:** `mstate` / `prodlim` (multi-state transition probabilities and
> state-occupancy curves); `survival` recurrent-event tooling
> (`survfit(Surv(t0, t1, status) ~ 1)` mean cumulative function / Andersen–Gill
> counting-process form); closed-form on an analytic illness-death DGP.

## Goal

Extend survatr beyond the **single terminal event** assumption to two related
event structures:

1. **Recurrent events** — an individual can experience the same event type
   multiple times (hospitalizations, relapses). Estimand: the counterfactual
   **mean cumulative function** (MCF) `μ^a(t) = E[N^a(t)]`, the expected number
   of events by `t`, built from a per-period **recurrence hazard** that allows
   the at-risk indicator to *reset* after each event rather than being absorbed.
2. **Multi-state / illness-death models** — individuals move among a finite set
   of states via **transition-specific hazards**; the estimands are
   counterfactual **state-occupancy probabilities** `P^a(state = s at t)` and
   transition probabilities, the multi-state generalization of the chunk-7 CIF
   (which is itself the special case of an absorbing-competing-risks multi-state
   model).

This is the **largest lift** in the roadmap. The package's rectangular,
single-terminal-event person-period assumption (the `survatr_ragged_pp` +
absorbing-event risk-set contract) **breaks** here: the risk set no longer
drops permanently at first event, and "the outcome" is no longer one indicator
per row. This chunk is as much a data-model redesign as an estimator addition.
We are honest about that up front and lay out the data-model changes and the
open design questions explicitly rather than pretending it slots into the
existing spine.

## The math

### Recurrent events — mean cumulative function

Per-period recurrence hazard (intensity) `ρ(t | A, L)` on rows where the
individual is **at risk for a (possibly repeat) event**, fit by pooled logistic
on a person-period table whose risk-set indicator does **not** absorb at first
event:

```
logit ρ(t | A, L) = alpha(t) + beta_A A + beta_L L           [at-risk, recurrence-resetting]
```

The counterfactual MCF accumulates the expected per-period event count
(optionally weighted by survival to `t` when a terminal event is also present —
the Cook–Lawless distinction between the marginal MCF and the
terminal-conditional MCF, an open design question below):

```
μ^a_i(t) = Σ_{k <= t}  ρ̂^a_{i,k}          (no-terminal-event case)
μ^a(t)   = (1/n) Σ_i μ^a_i(t)
```

The cumulative **sum** of intensities replaces the cumulative **product** of
survival-hazards. Variance is the chunk-3 delta chain, but the per-individual
functional is the running sum (the gradient is `∂μ/∂ρ_k = 1`) rather than the
running product — a simpler delta than survival, on the same IF spine.

### Multi-state — transition-specific hazards + state occupancy

For a state space `S = {1, …, R}` with allowed transitions `r → r'`, fit a
**transition-specific hazard** `h^{(r→r')}(t | A, L)` on rows where the
individual currently occupies state `r` (a cause-specific hazard generalized
from "leave alive" to "leave state `r` for state `r'`"):

```
logit h^{(r→r')}(t | A, L) = alpha_{r,r'}(t) + beta_{A,r,r'} A + beta_{L,r,r'} L
```

Counterfactual state-occupancy probabilities follow the discrete-time
**Aalen–Johansen** product-integral: with the per-period transition-probability
matrix `P^a(k)` (off-diagonals = transition hazards, diagonal = stay
probability), the occupancy vector is the running matrix product

```
π^a(t) = π^a(0) · ∏_{k <= t} P^a(k)
```

The chunk-7 competing-risks CIF is exactly this with one transient state and `J`
absorbing states — so the chunk-7 stacked-EE cause-specific machinery
**generalizes**: the parallel cause-specific fits become parallel
transition-specific fits, and the CIF sum becomes the Aalen–Johansen product.
Variance stacks the transition-specific scores; the per-individual IF is again
`n × |t-grid|` (now per occupancy probability / per state), full time-covariance
`crossprod(IF) / n²`.

## Data-model changes required (the honest part)

The single-terminal-event PP contract does not survive. Concretely:

- **Risk set no longer absorbs at first event.** `risk_set.R`'s
  `.survatr_prev_event` cumsum (drop rows at/after first event) is wrong for
  recurrent events; the at-risk indicator must **reset** after each event. A new
  recurrence-resetting risk-set builder is needed.
- **"Outcome" is no longer one 0/1 indicator per row.** For multi-state, each
  at-risk row needs a **from-state** and (on a transition) a **to-state**; for
  recurrent events, the row carries an event count / repeat-event indicator.
  `prepare_data.R`'s `survatr_bad_indicator` (0/1 only) and the rectangular
  guard must be replaced by a state-table validator.
- **A state / transition specification** is required from the user — the state
  space, the allowed transitions, and the initial-state distribution. This is a
  genuinely new input object, not a tweak to `confounders =` / `outcome =`.
- **The estimand shape gains a `state` (or `transition`) dimension** on
  `survatr_result`, analogous to chunk-7's `cause` dimension but richer (a
  vector of occupancy probabilities per time, not a single CIF).

These changes are why this is a separate large chunk and not a fold-in.

## Deliverables

### New R files
- `R/recurrent_multistate.R` — the recurrence-resetting risk-set builder; the
  transition-specific parallel fits; the MCF (cumulative-sum) estimand; the
  Aalen–Johansen state-occupancy product; the stacked-EE sandwich generalizing
  chunk-7; the state-table validator and the state/transition specification
  object.

### Updated R files
- `R/surv_fit.R` — accept a `states =` / `transitions =` specification (and a
  recurrent-event mode flag); route to the new risk-set builder; set a new
  track / event-structure tag on the `survatr_fit`.
- `R/risk_set.R` — add the recurrence-resetting at-risk logic (do not overwrite
  the absorbing-event path used by chunks 1–7; branch on the event structure).
- `R/prepare_data.R` — state-table validator replacing the 0/1-indicator +
  rectangular guards for the multi-state / recurrent path; classed aborts
  (`survatr_bad_state_table`, `survatr_bad_transition_spec`).
- `R/contrast.R` — MCF / state-occupancy / transition-probability estimands;
  the `state` / `transition` result dimension.
- `R/variance_if_survival.R` — cumulative-sum delta chain (MCF) and the
  matrix-product (Aalen–Johansen) delta chain (state occupancy).
- `R/plot.R` / `R/tidy.R` / `R/print.R` — render stacked state-occupancy curves
  and the MCF.

### Tests (`tests/testthat/`)
- `test-recurrent-events.R` — analytic recurrent-event DGP: the counterfactual
  MCF matches the closed-form `Σ_k ρ^a(k)`; vs `survival` MCF
  (`survfit(Surv(t0, t1, status) ~ 1)`) on the unadjusted curve; the
  recurrence-resetting risk set is asserted (the at-risk indicator returns to 1
  after an event).
- `test-multistate.R` — analytic illness-death DGP: state-occupancy
  probabilities match the closed-form Aalen–Johansen product; the rows sum to 1
  across states at each `t` (the occupancy-probability identity); the chunk-7
  CIF is recovered exactly as the one-transient-two-absorbing special case.
- `test-multistate-sandwich.R` — coverage simulation; transition-probability
  contrast CI covers 0 under a no-effect DGP; cross-validate the occupancy
  curve against `mstate` / `prodlim`.
- `helper-multistate-oracle.R` — `mstate` / `prodlim` comparators + the analytic
  Aalen–Johansen reference.

## API contract

```r
# Recurrent events — counterfactual mean cumulative function
fit <- surv_fit(
  data, outcome = "event",          # event count / repeat-event indicator
  treatment = "A", confounders = ~ L,
  id = "id", time = "t",
  recurrent = TRUE,                 # NEW — recurrence-resetting risk set
  estimator = "gcomp"
)
result <- contrast(fit, interventions = list(a1 = causatr::static(1),
                                             a0 = causatr::static(0)),
                   times = ..., type = "mcf_difference")
# result$estimates: intervention | time | mcf_hat | se | ci_*

# Multi-state — state occupancy via transition-specific hazards
fit <- surv_fit(
  data, state = "state",            # NEW — current-state column
  treatment = "A", confounders = ~ L, id = "id", time = "t",
  transitions = list(c(1, 2), c(1, 3), c(2, 3)),   # NEW — allowed r -> r'
  estimator = "gcomp"
)
result <- contrast(fit, ..., type = "state_occupancy", state = 2)
# result$estimates: intervention | state | time | prob_hat | se | ci_*
```

## Behaviour rules (non-negotiable — see hard-rules.md)

- **Recurrent events use a recurrence-RESETTING risk set** — the at-risk
  indicator returns to 1 after each event; the absorbing-event
  `.survatr_prev_event` cumsum (chunks 1–7) is wrong here and must not be
  reused. Branch on the event structure; do not break the absorbing path.
- **The MCF is a cumulative SUM of intensities**, not a cumulative product of
  survival-hazards (recurrent events are counts, not a single survival).
- **Multi-state occupancy uses the Aalen–Johansen product-integral** of
  per-period transition-probability matrices; the chunk-7 CIF is the
  one-transient-state special case and must be recovered exactly.
- **State-occupancy probabilities sum to 1 across states at each `t`** (the
  occupancy identity), the multi-state analogue of the chunk-7
  `Σ_j F^(j)(t) + S(t) = 1` identity — must hold numerically.
- **The transition-specific stacked-EE sandwich generalizes chunk 7** — reuse
  the parallel-fit + stacked-score structure; do not reimplement.
- **Fine–Gray and any subdistribution / pseudo-value approach remain out of
  scope** (the chunk-7 rule carries forward to this richer setting).
- **The single-terminal-event guards do not apply on this path** — the
  rectangular + 0/1-indicator validators are replaced by the state-table
  validator; the absorbing-event path keeps its guards unchanged.

## Open design questions (flag explicitly — resolve at implementation time)

- **Recurrent + terminal event.** Marginal MCF (Cook–Lawless) vs
  terminal-conditional MCF (weight intensities by survival to `t`). Which is the
  default estimand, and is the choice user-facing?
- **Within-id event correlation.** Recurrent-event variance almost always needs
  the chunk-13 cluster-robust SE keyed on id (events within id are correlated).
  Is clustered SE the default for recurrent events, or opt-in?
- **State specification ergonomics.** A transition list vs a transition matrix
  vs an `mstate`-style `transMat`; how much to lean on an existing convention
  vs a survatr-native object.
- **Reversible transitions** (recovery, `2 → 1`). Does v1 of this chunk allow
  reversible illness-death, or restrict to progressive (no back-transitions)
  multi-state first?
- **Time-varying treatment in multi-state** (Track B × multi-state). Almost
  certainly a separate future chunk; confirm it is deferred, not folded in.
- **Identifiability / positivity per transition.** Sparse late-time transition
  risk sets are worse than in single-event survival; how aggressively should
  `diagnose()` (chunk 10) flag empty transition cells?

## Non-goals (deferred)
- **Semi-competing risks** (illness-death where the non-terminal event is
  subject to the terminal event) as a distinct estimand class — adjacent, but
  its own modelling story; defer.
- **Frailty / random-effect recurrence models** — out of scope; recurrence
  correlation is handled by clustered SE (chunk 13), not a frailty term.
- **Track B (time-varying treatment) × multi-state** — a separate future lift.
- **Continuous-time multi-state** (`msm`-style transition-intensity matrices on
  continuous time) — survatr is discrete-time pooled-logistic by contract.

## Dependencies & composition
- Chunks 2, 3, 7 (the CIF / cause-specific machinery is the direct ancestor of
  transition-specific hazards). Chunk 13 (cluster-robust SE) is effectively a
  companion for recurrent events. No new causatr internals beyond
  `causatr:::prepare_model_if()` per transition-specific model (already imported
  by chunks 3/7); confirm exact symbol names against the current causatr `R/` at
  implementation time and pin them in `test-causatr-integration.R`.

## Acceptance checklist
- [ ] Recurrence-resetting risk set built (at-risk returns to 1 after an event);
      counterfactual MCF matches the analytic `Σ_k ρ^a(k)` and `survival`'s MCF.
- [ ] Multi-state state-occupancy matches the analytic Aalen–Johansen product;
      occupancy probabilities sum to 1 across states at each `t`.
- [ ] The chunk-7 CIF is recovered exactly as the one-transient-state special
      case.
- [ ] Transition-specific stacked-EE sandwich covers nominally; cross-validated
      against `mstate` / `prodlim`.
- [ ] State-table / transition-spec validators abort on malformed input; the
      absorbing-event path retains its guards unchanged.
- [ ] Open design questions resolved + recorded in CLAUDE.md scope.
- [ ] `FEATURE_COVERAGE_MATRIX.md` + handoff §10 + CLAUDE.md updated.
