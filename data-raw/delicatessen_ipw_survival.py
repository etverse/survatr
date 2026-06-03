"""
Reference values for the IPW survival weighted-hazard MSM via delicatessen
stacked M-estimation (Zivich 2022). Independent analytic-sandwich oracle for
survatr's `estimator = "ipw"` + `ci_method = "sandwich"` path.

Binary point treatment, discrete-time survival on a rectangular person-period
grid. Estimand: counterfactual survival S^a(t) under static(a), a in {0, 1},
with stabilized observed-treatment weights w_i = f(A_i) / f(A_i | L_i).

Stacked estimating-equation system (one observation = one id):
  1. Propensity (alpha), on the baseline row per id:
       Z_i' (A_i - expit(Z_i alpha)),  Z_i = [1, L_i]
  2. Weighted pooled-logistic hazard MSM (beta), summed over each id's
     at-risk person-period rows:
       sum_k atrisk_ik * X_ik' * w_i(alpha) * (Y_ik - expit(X_ik beta))
     X_ik = [1, factor(t) dummies, A_i]  (the survatr `~ factor(t) + A` design).
  3. Survival functional (theta_{a,t}), per id, over ALL periods k <= t:
       prod_{k<=t} (1 - expit(X_ik^{A=a} beta)) - theta_{a,t}

The stabilizing numerator f(A_i) is the FIXED marginal mean of A (pm), matching
survatr: only the denominator propensity is propagated into the variance (the
marginal model's estimation is omitted, a conservative choice -- Robins 1999).
This is what makes the delicatessen sandwich comparable to survatr's stacked
sandwich rather than narrower.

Usage (delicatessen lives in causatr's data-raw venv):
    ../causatr/data-raw/zepid_venv/bin/python delicatessen_ipw_survival.py

Reads  tests/testthat/fixtures/python/ipw_survival_data.csv  (shared with R)
Writes tests/testthat/fixtures/python/ipw_survival_delicatessen.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator

TIMES = [2, 3, 4, 5]  # survival horizons to pin


def build_design(periods, K, a):
    """One person-period design row [1, factor(t)_2..K dummies, a] per period."""
    n = len(periods)
    X = np.zeros((n, K + 1))
    X[:, 0] = 1.0
    for k in range(2, K + 1):
        X[:, k - 1] = (periods == k).astype(float)
    X[:, K] = a
    return X


def ipw_survival_reference(df):
    ids = np.sort(df["id"].unique())
    n = len(ids)
    K = int(df["t"].max())
    p_beta = K + 1  # intercept + (K-1) time dummies + A

    # Per-id arrays. L and A are baseline-constant within id.
    L = np.array([df.loc[df["id"] == i, "L"].iloc[0] for i in ids])
    A = np.array([df.loc[df["id"] == i, "A"].iloc[0] for i in ids], dtype=float)
    # Full K-period stacks (rectangular grid), ordered by period.
    Y = np.zeros((n, K))
    atrisk = np.zeros((n, K))
    for r, i in enumerate(ids):
        sub = df[df["id"] == i].sort_values("t")
        Y[r, :] = sub["Y"].values
        atrisk[r, :] = sub["at_risk"].values
    periods = np.arange(1, K + 1)

    # Observed-A design (for the hazard score) and counterfactual designs.
    Xobs = np.stack([build_design(periods, K, A[r]) for r in range(n)])  # n,K,p
    X_a = {a: build_design(periods, K, float(a)) for a in (0, 1)}  # K x p each

    Z = np.column_stack([np.ones(n), L])  # baseline propensity design
    p_alpha = Z.shape[1]
    pm = A.mean()  # FIXED stabilizing numerator P(A = 1)
    f_num = np.where(A == 1, pm, 1 - pm)

    n_surv = 2 * len(TIMES)
    # theta = [alpha(p_alpha), beta(p_beta), S1_t..., S0_t...]
    i_beta = p_alpha
    i_surv = p_alpha + p_beta

    def psi(theta):
        alpha = theta[:p_alpha]
        beta = theta[i_beta:i_beta + p_beta]
        surv = theta[i_surv:]

        # --- propensity (per id) ---
        pi_val = expit(Z @ alpha)
        ee_alpha = Z.T * (A - pi_val)  # p_alpha x n

        # --- stabilized weight (denominator depends on alpha) ---
        f_den = np.where(A == 1, pi_val, 1 - pi_val)
        w = f_num / f_den  # n

        # --- weighted hazard MSM score, summed over at-risk PP rows per id ---
        ee_beta = np.zeros((p_beta, n))
        h_obs = expit(np.einsum("nkp,p->nk", Xobs, beta))  # n x K
        resid = (Y - h_obs) * atrisk  # zero off the risk set
        # sum_k X_ik * resid_ik, scaled by w_i
        ee_beta = np.einsum("nkp,nk->pn", Xobs, resid) * w[np.newaxis, :]

        # --- survival functional per id, per (a, t) ---
        ee_surv = np.zeros((n_surv, n))
        row = 0
        for a in (1, 0):
            h_a = expit(X_a[a] @ beta)  # length K (same for every id; A const)
            surv_cum = np.cumprod(1 - h_a)  # S^a(t) for t = 1..K (id-invariant)
            for t in TIMES:
                # S^a_i(t) is identical across ids here (marginal MSM, no EM),
                # but written per id so a future effect-modifier reuses it.
                ee_surv[row, :] = surv_cum[t - 1] - surv[row]
                row += 1

        return np.vstack([ee_alpha, ee_beta, ee_surv])

    # --- initial values ---
    def nll(b, Xd, yd, wd):
        p = np.clip(expit(Xd @ b), 1e-10, 1 - 1e-10)
        return -np.sum(wd * (yd * np.log(p) + (1 - yd) * np.log(1 - p)))

    a_init = minimize(
        lambda b: nll(b, Z, A, np.ones(n)), np.zeros(p_alpha), method="BFGS"
    ).x
    pi0 = expit(Z @ a_init)
    w0 = f_num / np.where(A == 1, pi0, 1 - pi0)
    # at-risk rows flattened for the beta init fit
    Xflat = Xobs.reshape(n * K, p_beta)
    yflat = Y.reshape(n * K)
    wflat = (atrisk * w0[:, None]).reshape(n * K)
    b_init = minimize(
        lambda b: nll(b, Xflat, yflat, wflat), np.zeros(p_beta), method="BFGS"
    ).x
    surv_init = []
    for a in (1, 0):
        sc = np.cumprod(1 - expit(X_a[a] @ b_init))
        surv_init += [sc[t - 1] for t in TIMES]
    init = np.concatenate([a_init, b_init, surv_init])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")
    theta = mest.theta
    vcov = mest.variance

    surv_idx = {}
    row = 0
    for a in (1, 0):
        for t in TIMES:
            surv_idx[(a, t)] = i_surv + row
            row += 1

    out = []
    for a in (1, 0):
        for t in TIMES:
            j = surv_idx[(a, t)]
            out.append(
                dict(quantity=f"S{a}", time=t,
                     estimate=theta[j], se=np.sqrt(vcov[j, j]))
            )
    # RD(t) = S1(t) - S0(t) (survatr "a0 vs a1" = risk0 - risk1 = S1 - S0)
    for t in TIMES:
        j1, j0 = surv_idx[(1, t)], surv_idx[(0, t)]
        est = theta[j1] - theta[j0]
        var = vcov[j1, j1] + vcov[j0, j0] - 2 * vcov[j1, j0]
        out.append(dict(quantity="RD", time=t, estimate=est, se=np.sqrt(var)))
    return pd.DataFrame(out)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fx = os.path.normpath(
        os.path.join(here, "..", "tests", "testthat", "fixtures", "python")
    )
    df = pd.read_csv(os.path.join(fx, "ipw_survival_data.csv"))
    print(f"Loaded {len(df)} rows, {df['id'].nunique()} ids, K={df['t'].max()}")
    res = ipw_survival_reference(df)
    res["estimate"] = res["estimate"].round(8)
    res["se"] = res["se"].round(8)
    print(res.to_string(index=False))
    res.to_csv(os.path.join(fx, "ipw_survival_delicatessen.csv"), index=False)
    print("wrote ipw_survival_delicatessen.csv")
