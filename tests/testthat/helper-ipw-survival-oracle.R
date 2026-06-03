## Reusable lmtp comparator for the IPW survival point estimate.
##
## `lmtp::lmtp_tmle(outcome_type = "survival")` is a TMLE estimator of the
## counterfactual survival functional S^a(tau) under a static intervention; with
## a correctly specified treatment / outcome it is consistent and serves as an
## external point-estimate oracle for survatr's IPW weighted-MSM curve. (Its
## EIF-based SE is NOT comparable to the M-estimation sandwich, so this is a
## point oracle only -- the sandwich is validated against the bootstrap.)
##
## lmtp's survival estimand is the survival probability at the LAST supplied
## outcome period, so the full curve is recovered by calling it once per horizon
## with the outcome / censoring vectors truncated to that horizon. The scalar
## estimate lives in the S7 `ife` object at `@theta`.
##
## Reshapes a rectangular person-period `data.table` (columns id, t, A, L, Y)
## into the wide one-row-per-id layout, then for each horizon in `times` fits
## the static(value) intervention adjusting for the baseline confounder L.
## Returns the lmtp S(t) vector (one entry per requested horizon), or NULL if
## lmtp is unavailable or any horizon fails to fit.
lmtp_ipw_survival_oracle <- function(
  dt,
  value,
  confounder = "L",
  times = NULL
) {
  if (!requireNamespace("lmtp", quietly = TRUE)) {
    return(NULL)
  }
  dt <- data.table::as.data.table(dt)
  K <- length(unique(dt$t))
  if (is.null(times)) {
    times <- seq_len(K)
  }

  ## Wide reshape: one row per id with A, the baseline confounder(s), and the
  ## per-period event indicators Y_1..Y_K. `...` keeps id-constant columns (L).
  wide <- data.table::dcast(dt, id + A + ... ~ t, value.var = "Y")
  data.table::setnames(wide, as.character(seq_len(K)), paste0("Y_", seq_len(K)))
  ## lmtp's survival outcome is the cumulative event indicator.
  for (k in 2:K) {
    wide[[paste0("Y_", k)]] <- pmax(
      wide[[paste0("Y_", k - 1L)]],
      wide[[paste0("Y_", k)]]
    )
  }
  ## Always-observed censoring indicators (the DGP has no censoring).
  for (k in seq_len(K)) {
    wide[[paste0("C_", k)]] <- 1L
  }
  wide_df <- as.data.frame(wide)

  out <- vapply(
    times,
    function(tau) {
      fit_lmtp <- tryCatch(
        lmtp::lmtp_tmle(
          data = wide_df,
          trt = "A",
          outcome = paste0("Y_", seq_len(tau)),
          baseline = confounder,
          cens = paste0("C_", seq_len(tau)),
          shift = function(data, trt) rep(value, nrow(data)),
          outcome_type = "survival",
          folds = 1L
        ),
        error = function(e) NULL
      )
      if (is.null(fit_lmtp)) {
        return(NA_real_)
      }
      ## S7 `ife` object: the survival point estimate at horizon `tau` is `@x`.
      as.numeric(fit_lmtp$estimate@x)
    },
    numeric(1L)
  )

  if (anyNA(out)) {
    return(NULL)
  }
  out
}
