# surv_fit rejects matching estimator with pointer error

    Code
      surv_fit(dt, "Y", "A", ~L, "id", "t", estimator = "matching")
    Condition
      Error in `surv_fit()`:
      ! Matching + survival is out of scope for survatr.
      i Use `survival::coxph(..., weights = match_weights, cluster = subclass)` directly on the `MatchIt` output.

# print.survatr_fit emits a stable banner

    Code
      print(fit)
    Output
      <survatr_fit>
        Track:       A
        Estimator:   gcomp
        Family:      binomial
        Outcome:     Y
        Treatment:   A
        ID:          id
        Time:        t
        Censoring:   cens
        N:           5 individuals, 15 PP rows (10 at risk)
        Time grid:   [1, 3] (3 unique times)

