# prepare_pp_data rejects missing columns

    Code
      prepare_pp_data(dt, "Y", "A", "id", "t", censoring = "cens_missing")
    Condition
      Error:
      ! Column(s) `cens_missing` not found in `data`.

# prepare_pp_data rejects wide (one row per id) input

    Code
      prepare_pp_data(wide, "Y", "A", "id", "t")
    Condition
      Error:
      ! Data does not appear to be in person-period format: 5 of 5 unique `id` values have only a single row. Use `causatr::to_person_period()` to convert wide data to long before calling `surv_fit()`.

