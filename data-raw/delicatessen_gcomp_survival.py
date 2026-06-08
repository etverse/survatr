"""
Reference values for the Track A g-computation (gcomp) survival sandwich via
delicatessen stacked M-estimation (Zivich 2022). Independent analytic-sandwich
oracle for survatr's `estimator = "gcomp"` + `ci_method = "sandwich"` path --
the foundational variance chunk that previously had only a closed-form analytic
SE + empirical-SD + bootstrap check (no independent M-estimator).

Binary point treatment, discrete-time survival on a rectangular person-period
grid. g-computation (standardization over the empirical L distribution):
S^a(t) = mean_i prod_{k<=t} (1 - expit(X_ik^{A=a} beta)),
X_ik = [1, factor(t) dummies, A_i, L_i] (the survatr ~ factor(t) + A + L design).
Confounding is handled by the outcome model (L in the design), NOT by weights.

Stacked estimating-equation system (one observation = one id):
  1. Pooled-logistic hazard (beta), summed over each id's at-risk PP rows:
       sum_k atrisk_ik * X_ik' * (Y_ik - expit(X_ik beta))   (plain binomial)
  2. Survival functional (theta_{a,t}), per id, over ALL periods k <= t:
       prod_{k<=t} (1 - expit(X_ik^{A=a} beta)) - theta_{a,t}

Risk difference, risk ratio, and RMST (and their SEs) are derived post hoc from
the survival-theta point estimates + covariance block (delta method), exactly as
survatr aggregates its survival influence function -- so the comparison pins the
whole gcomp estimand surface, not just S^a(t). The shared person-period fixture
is the same one the IPW oracle reads (the gcomp and IPW estimators target the
same counterfactual S^a(t) on the same confounded data).

Usage (delicatessen lives in causatr's data-raw venv):
    ../causatr/data-raw/zepid_venv/bin/python delicatessen_gcomp_survival.py

Reads  tests/testthat/fixtures/python/ipw_survival_data.csv  (shared with R)
Writes tests/testthat/fixtures/python/gcomp_survival_delicatessen.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator

TIMES = [2, 3, 4, 5]  # survival horizons to pin (also the RMST/RD/RR grid)


def build_design(periods, K, a, L):
    """One PP design row [1, factor(t)_2..K, a, L] per period (A, L baseline)."""
    n = len(periods)
    X = np.zeros((n, K + 2))
    X[:, 0] = 1.0
    for k in range(2, K + 1):
        X[:, k - 1] = (periods == k).astype(float)
    X[:, K] = a
    X[:, K + 1] = L
    return X


def rmst_weights(times):
    """Trapezoidal Jacobian d RMST(t_j) / d [S(t_1..t_K)] (survatr rmst_weights).

    With dt[i] = t_i - t_{i-1}, t_0 = 0:
      W[j, i]  = dt[i]/2                    (i == 1)
               + dt[i+1]/2 if i < j         (interior)
      W[j, j] += dt[j]/2                     (last half-interval)
    The S(0)=1 constant dt[1]/2 is added separately (drops out of the delta).
    """
    K = len(times)
    dt = np.diff([0] + list(times))
    W = np.zeros((K, K))
    for j in range(K):
        for i in range(j + 1):
            W[j, i] += dt[i] / 2.0
            if i < j:
                W[j, i] += dt[i + 1] / 2.0
    return W, dt


def gcomp_reference(df):
    ids = np.sort(df["id"].unique())
    n = len(ids)
    K = int(df["t"].max())
    p_beta = K + 2  # intercept + (K-1) time dummies + A + L
    periods = np.arange(1, K + 1)

    L = np.array([df.loc[df["id"] == i, "L"].iloc[0] for i in ids])
    A = np.array([df.loc[df["id"] == i, "A"].iloc[0] for i in ids], dtype=float)
    Y = np.zeros((n, K))
    atrisk = np.zeros((n, K))
    for r, i in enumerate(ids):
        sub = df[df["id"] == i].sort_values("t")
        Y[r, :] = sub["Y"].values
        atrisk[r, :] = sub["at_risk"].values

    Xobs = np.stack([build_design(periods, K, A[r], L[r]) for r in range(n)])
    X_a = {
        a: np.stack([build_design(periods, K, float(a), L[r]) for r in range(n)])
        for a in (0, 1)
    }

    n_surv = 2 * len(TIMES)  # (a in {1,0}) x times
    i_surv = p_beta

    surv_idx = {}
    row = 0
    for a in (1, 0):
        for t in TIMES:
            surv_idx[(a, t)] = i_surv + row
            row += 1

    def psi(theta):
        beta = theta[:p_beta]
        surv = theta[i_surv:]

        # --- pooled-logistic hazard score, summed over at-risk PP rows ---
        h_obs = expit(np.einsum("nkp,p->nk", Xobs, beta))
        resid = (Y - h_obs) * atrisk
        ee_beta = np.einsum("nkp,nk->pn", Xobs, resid)

        # --- per-id standardized survival functional ---
        ee_surv = np.zeros((n_surv, n))
        for a in (1, 0):
            h_a = expit(np.einsum("nkp,p->nk", X_a[a], beta))  # n x K
            surv_cum = np.cumprod(1 - h_a, axis=1)  # S^a_i(t), t = 1..K
            for t in TIMES:
                k = surv_idx[(a, t)] - i_surv
                ee_surv[k, :] = surv_cum[:, t - 1] - surv[k]

        return np.vstack([ee_beta, ee_surv])

    # --- initial values: plain logistic hazard fit on at-risk rows ---
    def nll(b):
        Xf = Xobs.reshape(n * K, p_beta)
        yf = Y.reshape(n * K)
        wf = atrisk.reshape(n * K)
        pr = np.clip(expit(Xf @ b), 1e-10, 1 - 1e-10)
        return -np.sum(wf * (yf * np.log(pr) + (1 - yf) * np.log(1 - pr)))

    b_init = minimize(nll, np.zeros(p_beta), method="BFGS").x
    surv_init = []
    for a in (1, 0):
        sc = np.cumprod(1 - expit(np.einsum("nkp,p->nk", X_a[a], b_init)), axis=1)
        surv_init += [sc[:, t - 1].mean() for t in TIMES]
    init = np.concatenate([b_init, surv_init])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")
    theta = mest.theta
    vcov = mest.variance

    def sd(j):
        return float(np.sqrt(vcov[j, j]))

    out = []
    for a in (1, 0):
        for t in TIMES:
            j = surv_idx[(a, t)]
            out.append(dict(quantity=f"S{a}", time=t,
                            estimate=theta[j], se=sd(j)))
    # RD(t) = S1 - S0 (survatr "a0 vs a1" risk-difference = risk0 - risk1 = S1 - S0)
    for t in TIMES:
        j1, j0 = surv_idx[(1, t)], surv_idx[(0, t)]
        est = theta[j1] - theta[j0]
        var = vcov[j1, j1] + vcov[j0, j0] - 2 * vcov[j1, j0]
        out.append(dict(quantity="RD", time=t, estimate=est, se=np.sqrt(var)))
    # RR(t) on risk = 1 - S, survatr "a0 vs a1" = risk_a0/risk_a1 = r0/r1;
    # delta method on the S covariance.
    for t in TIMES:
        j1, j0 = surv_idx[(1, t)], surv_idx[(0, t)]
        r1, r0 = 1 - theta[j1], 1 - theta[j0]
        rr = r0 / r1
        g = np.zeros(len(theta))
        g[j1] = r0 / r1**2     # d RR / d S1
        g[j0] = -1.0 / r1      # d RR / d S0
        se = float(np.sqrt(g @ vcov @ g))
        out.append(dict(quantity="RR", time=t, estimate=rr, se=se))
    # RMST^a(t_j) = W_j . [S^a(t_1..t_K)] + dt[1]/2; SE from the S covariance.
    W, dt = rmst_weights(TIMES)
    for a in (1, 0):
        js = [surv_idx[(a, t)] for t in TIMES]
        s_vec = theta[js]
        cov_ss = vcov[np.ix_(js, js)]
        for jrow, t in enumerate(TIMES):
            w = W[jrow]
            rmst = float(w @ s_vec + dt[0] / 2.0)
            se = float(np.sqrt(w @ cov_ss @ w))
            out.append(dict(quantity=f"RMST{a}", time=t,
                            estimate=rmst, se=se))
    # RMST difference, survatr "a0 vs a1" = rmst_a0 - rmst_a1.
    js1 = [surv_idx[(1, t)] for t in TIMES]
    js0 = [surv_idx[(0, t)] for t in TIMES]
    for jrow, t in enumerate(TIMES):
        w = W[jrow]
        g = np.zeros(len(theta))
        g[js0] += w
        g[js1] -= w
        est = float(g @ theta)
        se = float(np.sqrt(g @ vcov @ g))
        out.append(dict(quantity="RMSTdiff", time=t, estimate=est, se=se))
    return pd.DataFrame(out)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fx = os.path.normpath(
        os.path.join(here, "..", "tests", "testthat", "fixtures", "python")
    )
    df = pd.read_csv(os.path.join(fx, "ipw_survival_data.csv"))
    print(f"Loaded {len(df)} rows, {df['id'].nunique()} ids, K={df['t'].max()}")
    res = gcomp_reference(df)
    res["estimate"] = res["estimate"].round(8)
    res["se"] = res["se"].round(8)
    print(res.to_string(index=False))
    res.to_csv(os.path.join(fx, "gcomp_survival_delicatessen.csv"), index=False)
    print("wrote gcomp_survival_delicatessen.csv")
