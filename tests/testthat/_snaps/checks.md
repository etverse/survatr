# check_weights rejects non-numeric / mis-sized / NA / Inf / negative

    Code
      check_weights(c("a", "b", "c"), 3L)
    Condition
      Error:
      ! `weights` must be numeric.

---

    Code
      check_weights(c(1, 2, 3), 5L)
    Condition
      Error:
      ! `weights` must have length equal to `nrow(data)` (5), got 3.

---

    Code
      check_weights(c(1, NA, 3), 3L)
    Condition
      Error:
      ! `weights` contains 1 missing value(s). Drop those rows or impute before calling `surv_fit()`.

---

    Code
      check_weights(c(1, Inf, 3), 3L)
    Condition
      Error:
      ! `weights` contains non-finite value(s) (Inf / NaN).

---

    Code
      check_weights(c(1, NaN, 3), 3L)
    Condition
      Error:
      ! `weights` contains 1 missing value(s). Drop those rows or impute before calling `surv_fit()`.

---

    Code
      check_weights(c(1, -2, 3), 3L)
    Condition
      Error:
      ! `weights` must be non-negative.

# check_dots_na_action rejects na.exclude (function and string)

    Code
      check_dots_na_action(na.action = stats::na.exclude)
    Condition
      Error:
      ! `na.action` must be `na.omit` (default) or `na.fail`.
      i survatr builds its own row-alignment bookkeeping from `fit_rows` and the fitted model's `na.action` attribute. `na.exclude` pads working residuals with NA and silently corrupts the sandwich variance.
      i Drop NA rows before calling `surv_fit()` or use `na.action = na.omit`.

# check_reserved_cols rejects .survatr_prev_event / .survatr_prev_cens

    Code
      check_reserved_cols(bad_both)
    Condition
      Error:
      ! Column name(s) `.survatr_prev_event`, `.survatr_prev_cens` are reserved by survatr internals. Rename the column(s) in your input data before calling `surv_fit()`.

