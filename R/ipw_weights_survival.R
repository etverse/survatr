## IPW survival weight helpers: stabilized density-ratio weight composition,
## the per-id -> person-period broadcast, and the alpha-closure used by the
## stacked-EE sandwich. Split out of `ipw_survival.R` (the orchestrator) to keep
## each file under the ~300-line soft cap. causatr supplies the density,
## positivity, and truncation primitives; survatr composes the stabilized ratio.

#' Compute stabilized density-ratio weights at the observed treatment
#'
#' @description
#' Form the per-id stabilized IPW weight `w_i = f(A_i) / f(A_i | L_i)` from a
#' fitted denominator (full) and numerator (marginal) treatment model, evaluated
#' at the observed treatment, and optionally winsorize at a fixed cutoff.
#'
#' @details
#' Both densities come from `causatr:::evaluate_density()`. The denominator is
#' guarded for positivity via causatr's `check_density_positivity()` (a zero
#' conditional density is a hard positivity violation). Winsorization reuses
#' `causatr:::truncate_weights()` for the point weights but the resolved
#' threshold is captured and returned so the sandwich's weight closure can clip
#' at the **same fixed** cutoff (re-quantiling inside the numerical-derivative
#' closure would make the weight non-smooth in `alpha`).
#'
#' @param full_tm The full `A ~ L` `causatr_treatment_model` (denominator).
#' @param marginal_model The marginal `A ~ 1` `causatr_treatment_model`
#'   (numerator).
#' @param baseline One-row-per-id `data.table` (baseline rows).
#' @param treatment Treatment column name.
#' @param trim Numeric scalar in `(0, 1]` or `NULL`. `NULL` / `>= 1` means no
#'   winsorization.
#'
#' @returns A list with `weights` (numeric, length = number of baseline fit
#'   rows) and `trim_threshold` (the fixed cutoff, or `NA_real_` when no
#'   truncation was applied).
#' @family survatr_fit functions
#' @noRd
compute_ipw_stabilized_weights <- function(
  full_tm,
  marginal_model,
  baseline,
  treatment,
  trim = NULL
) {
  ## Evaluate on the rows the treatment model was actually fit on. Under
  ## surv_fit()'s upfront NA rejection these are all baseline rows, but honour
  ## `fit_rows` so the densities and weights stay aligned with `X_prop`.
  fit_data <- baseline[full_tm$fit_rows]
  a_obs <- fit_data[[treatment]]

  f_den <- causatr:::evaluate_density(full_tm, a_obs, fit_data)
  causatr:::check_density_positivity(f_den, "IPW stabilized weights")
  f_num <- causatr:::evaluate_density(marginal_model, a_obs, fit_data)

  w <- f_num / f_den

  ## Resolve the fixed winsorization cutoff once, at the point estimate, so the
  ## sandwich clips at the same value. `truncate_weights()` is a no-op for
  ## trim >= 1, matching the NULL case.
  trim_threshold <- NA_real_
  if (!is.null(trim) && trim < 1) {
    trim_threshold <- stats::quantile(w, trim, names = FALSE)
    w <- causatr:::truncate_weights(w, trim)
  }

  list(weights = w, trim_threshold = trim_threshold)
}

#' Broadcast per-id weights onto person-period rows
#'
#' @description
#' Map a per-id weight vector onto every person-period row by matching ids. The
#' weight is constant within id (a baseline quantity), so this is a pure lookup.
#'
#' @param w_by_id Numeric weights, one per baseline id, in `baseline_ids` order.
#' @param baseline_ids Vector of ids in the same order as `w_by_id`.
#' @param pp_ids Vector of ids, one per person-period row (the full PP column).
#'
#' @returns Numeric vector of length `length(pp_ids)`: each row gets its id's
#'   weight.
#' @family survatr_fit functions
#' @noRd
broadcast_weights_to_pp <- function(w_by_id, baseline_ids, pp_ids) {
  idx <- match(pp_ids, baseline_ids)
  if (anyNA(idx)) {
    ## Every PP id must have a baseline row (rectangular PP guarantees this).
    ## A missing match means a structural mismatch between the treatment-model
    ## frame and the PP frame -- fail loudly rather than recycle a wrong weight.
    rlang::abort(
      "Internal error: a person-period id has no baseline weight.",
      class = "survatr_ipw_broadcast_failed"
    )
  }
  w_by_id[idx]
}

#' Build the alpha-closure for the stabilized weight (binary treatment)
#'
#' @description
#' Return a closure `alpha -> per-id weight` used by the sandwich to
#' numerically differentiate the weighted hazard-MSM score with respect to the
#' propensity coefficients. Only the **denominator** `f(A | L; alpha)` depends
#' on `alpha`; the marginal numerator is captured fixed (omitting numerator
#' estimation from the IF is conservative -- Robins 1999; Hernán et al. 2000).
#'
#' @details
#' For a binary treatment with a logit propensity, `f(A_i | L_i; alpha) = p_i`
#' when `A_i = 1` and `1 - p_i` when `A_i = 0`, with `p_i = plogis(X_prop_i
#' alpha)`. The closure rebuilds this from the stored design `X_prop` (aliased
#' columns already dropped, aligned with `alpha_hat`), so the perturbed `alpha`
#' matches the bread / score column order. Winsorization, when active, clips at
#' the **fixed** threshold so the weight stays a smooth function of `alpha` away
#' from the (measure-zero) kink.
#'
#' @param X_prop The propensity design matrix (`full_tm$X_prop`); columns
#'   aligned with `full_tm$alpha_hat`.
#' @param a_obs Observed treatment values on the propensity fit rows (0/1).
#' @param fixed_numerator Marginal density `f(A_i)` on the fit rows, captured
#'   fixed.
#' @param trim_threshold Fixed winsorization cutoff, or `NA_real_` for none.
#'
#' @returns A function of `alpha` (numeric, length `ncol(X_prop)`) returning the
#'   per-id weight vector.
#' @family variance
#' @noRd
make_survival_weight_closure <- function(
  X_prop,
  a_obs,
  fixed_numerator,
  trim_threshold
) {
  force(X_prop)
  force(a_obs)
  force(fixed_numerator)
  force(trim_threshold)
  has_trim <- is.finite(trim_threshold)
  function(alpha) {
    eta <- as.numeric(X_prop %*% alpha)
    p <- stats::plogis(eta)
    ## Bernoulli density at the observed treatment.
    f_den <- ifelse(a_obs == 1, p, 1 - p)
    w <- fixed_numerator / f_den
    if (has_trim) {
      w <- pmin(w, trim_threshold)
    }
    w
  }
}
