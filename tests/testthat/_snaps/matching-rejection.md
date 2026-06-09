# rejects estimator = 'matching' with classed error and redirect

    Code
      surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "matching")
    Condition
      Error in `surv_fit()`:
      ! Matching + survival is out of scope for survatr.
      i Use `survival::coxph(..., weights = match_weights, cluster = subclass)` directly on the `MatchIt` output.

# rejects estimator = 'match' (alias) with classed error

    Code
      surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "match")
    Condition
      Error in `surv_fit()`:
      ! Matching + survival is out of scope for survatr.
      i Use `survival::coxph(..., weights = match_weights, cluster = subclass)` directly on the `MatchIt` output.

# rejects method = 'matching' in ... (causatr-style mis-call)

    Code
      surv_fit(dt, "Y", "A", ~L, "id", "t", method = "matching")
    Condition
      Error in `surv_fit()`:
      ! Matching + survival is out of scope for survatr.
      i Use `survival::coxph(..., weights = match_weights, cluster = subclass)` directly on the `MatchIt` output.

# rejects a `matchit` object passed as data

    Code
      surv_fit(mo, "Y", "A", ~L, "id", "t")
    Condition
      Error in `surv_fit()`:
      ! Matching + survival is out of scope for survatr.
      i Use `survival::coxph(..., weights = match_weights, cluster = subclass)` directly on the `MatchIt` output.

