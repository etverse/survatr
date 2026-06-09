## data-raw/nhefs_surv.R
## Assemble NHEFS monthly person-period survival data for the Ch. 17 replication.
## Source: causaldata package (Huntington-Klein 2021), which redistributes the
## NHEFS extract used in Hernán & Robins "Causal Inference: What If" (2024).
## Run once (not part of the build) to regenerate data/nhefs_surv.rda.
##
## Provenance: the NHEFS data are publicly distributed by Hernán & Robins at
## https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/ and are
## redistributed by the causaldata package under the same terms.

## ── Load and prepare ────────────────────────────────────────────────────────

nhefs_wide <- data.table::as.data.table(causaldata::nhefs)

## Compute the event time in monthly intervals (1 = January 1983, 120 = December
## 1992). For deaths (death == 1) the time is the month-of-death within the
## follow-up window: (yrdth - 83) converts the 2-digit year to years-since-1983,
## and *12 + modth yields the calendar month index (range 1..120).
## Non-deaths (survivors + the ~63 individuals lost to follow-up) are
## administratively censored at month 120.
nhefs_wide[,
  survtime := data.table::fifelse(
    death == 1L,
    as.integer((yrdth - 83L) * 12L + modth),
    120L
  )
]

## ── Person-period expansion (rectangular) ───────────────────────────────────
## surv_fit() requires a rectangular PP grid: every individual must have a row
## at every time period 1..120. The risk-set builder drops post-event rows from
## the hazard GLM fit (via the lagged cumsum of the event column), while the
## survival-curve prediction path uses all rows.

K <- 120L

surv_t <- nhefs_wide$survtime
died <- nhefs_wide$death == 1L
n_ids <- nrow(nhefs_wide)

## One row per (id, t) pair for t = 1..120.
seqn_long <- rep(nhefs_wide$seqn, each = K)
time_long <- rep(seq_len(K), n_ids)

## event = 1 only at the death month for individuals who died; 0 everywhere else.
event_vec <- integer(n_ids * K)
for (i in seq_len(n_ids)) {
  if (died[i]) {
    event_vec[(i - 1L) * K + surv_t[i]] <- 1L
  }
}

long_dt <- data.table::data.table(
  seqn = seqn_long,
  time = time_long,
  event = event_vec
)

## The nine confounders used in Hernán & Robins Ch. 17 (Table 17.1):
## sex, race, age (+ quadratic), education (5-cat), smokeintensity (+ quadratic),
## smokeyrs (+ quadratic), exercise (3-cat), active (3-cat), wt71 (+ quadratic).
baseline_cols <- c(
  "seqn",
  "qsmk",
  "sex",
  "race",
  "age",
  "education",
  "smokeintensity",
  "smokeyrs",
  "exercise",
  "active",
  "wt71"
)
baseline <- nhefs_wide[, ..baseline_cols]

## Left-join baseline onto the long grid: each (seqn, time) row inherits the
## individual's baseline covariates.
nhefs_surv <- merge(long_dt, baseline, by = "seqn", all.x = TRUE)
data.table::setcolorder(
  nhefs_surv,
  c(
    "seqn",
    "time",
    "event",
    "qsmk",
    "sex",
    "race",
    "age",
    "education",
    "smokeintensity",
    "smokeyrs",
    "exercise",
    "active",
    "wt71"
  )
)
data.table::setkeyv(nhefs_surv, c("seqn", "time"))

## ── Validate ─────────────────────────────────────────────────────────────────
stopifnot(nrow(nhefs_surv) == n_ids * K)
stopifnot(sum(nhefs_surv$event) == sum(died))
stopifnot(all(nhefs_surv$time >= 1L & nhefs_surv$time <= K))
stopifnot(!anyNA(nhefs_surv))

## ── Save ─────────────────────────────────────────────────────────────────────
dir.create("data", showWarnings = FALSE)
save(nhefs_surv, file = "data/nhefs_surv.rda", compress = "xz")
cat(
  "Saved data/nhefs_surv.rda:",
  nrow(nhefs_surv),
  "rows,",
  data.table::uniqueN(nhefs_surv$seqn),
  "individuals\n"
)
