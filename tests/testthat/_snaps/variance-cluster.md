# cluster validation aborts fire with classed errors

    Code
      contrast(fit_v, ivs, 1:5, "risk_difference", ci_method = "sandwich", cluster = "site")
    Condition
      Error in `contrast()`:
      ! `cluster = "site"` is not constant within `id`: 1 id(s) span more than one cluster.
      i Each individual must belong to exactly one cluster.

---

    Code
      contrast(fit_na, ivs, 1:5, "risk_difference", ci_method = "sandwich", cluster = "site")
    Condition
      Error in `contrast()`:
      ! `cluster = "site"` has NA values; every row must carry a cluster label.

---

    Code
      contrast(fit_1, ivs, 1:5, "risk_difference", ci_method = "sandwich", cluster = "site")
    Condition
      Error in `contrast()`:
      ! `cluster = "site"` resolves to a single cluster; cluster-robust variance needs at least two clusters.

---

    Code
      contrast(fit_ok, ivs, 1:5, "risk_difference", ci_method = "sandwich", cluster = "nope")
    Condition
      Error in `contrast()`:
      ! `cluster = "nope"` is not a column in the fit's person-period data.

