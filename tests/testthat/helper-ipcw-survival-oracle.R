## Reusable lmtp comparator for IPCW-corrected survival.
##
## lmtp::lmtp_tmle(outcome_type = "survival", cens = ...) is the external
## point-estimate oracle for survatr's IPCW-weighted hazard-MSM. With a
## correctly specified treatment, outcome, and censoring model, TMLE is
## consistent and converges to the same marginal counterfactual estimand as the
## IPCW g-computation curve. (lmtp's EIF-based SE is NOT comparable to the
## M-estimation sandwich; this is a point oracle only.)
##
## Convention difference: survatr uses C = 1 for "censored at period k",
## lmtp uses C_k = 1 for "observed at period k" (the inverse). This wrapper
## converts.
##
## Reshapes the rectangular person-period data.table (columns id, t, A, L, Y,
## C) into lmtp's wide one-row-per-id format with Y_1..Y_K (cumulative),
## C_1..C_K (not-censored), and fits TMLE for each horizon in `times`. Returns
## the S(t) vector (one per requested horizon), or NULL if lmtp is unavailable
## or a fit fails.
lmtp_ipcw_survival_oracle <- function(
  dt,
  value,
  confounder = "L",
  cens_col = "C",
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

  ## Wide reshape: one row per id; A and baseline confounder(s) stay as-is.
  Y_cols <- paste0("Y_", seq_len(K))
  C_cols <- paste0("C_", seq_len(K))
  wide <- data.table::dcast(dt, id + A + ... ~ t, value.var = c("Y", cens_col))
  ## Rename the Y_<t> columns to Y_1..Y_K.
  old_y <- paste0("Y_", seq_len(K))
  data.table::setnames(wide, old_y, Y_cols)
  old_c <- paste0(cens_col, "_", seq_len(K))
  data.table::setnames(wide, old_c, C_cols)

  ## lmtp survival outcome is cumulative: Y_k = 1 iff event by period k.
  for (k in 2L:K) {
    wide[[Y_cols[k]]] <- pmax(wide[[Y_cols[k - 1L]]], wide[[Y_cols[k]]])
  }

  ## lmtp censoring indicator: C_k = 1 means "observed" (not censored).
  ## survatr: C = 1 means censored. Invert.
  for (k in seq_len(K)) {
    wide[[C_cols[k]]] <- as.integer(wide[[C_cols[k]]] == 0L)
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
      as.numeric(fit_lmtp$estimate@x)
    },
    numeric(1L)
  )

  if (anyNA(out)) {
    return(NULL)
  }
  out
}
