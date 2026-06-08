"""
Reference values for competing-risks cause-specific hazards + CIF via
delicatessen stacked M-estimation (Zivich 2022). Independent analytic-sandwich
oracle for survatr's `competing =` + `ci_method = "sandwich"` path.

Two competing causes, binary point treatment, discrete-time on a rectangular
person-period grid. g-computation (standardization over the empirical L
distribution); no IPW weights. Estimand: cause-j cumulative incidence under
static(a), F^(j),a(t) = mean_i sum_{k<=t} S^a_i(k-1) h^(j),a_ik, with
S^a_i(k-1) = prod_{m<k} (1 - sum_j' h^(j'),a_im) and
h^(j),a_ik = expit(X_ik^{A=a} beta_j).

Stacked estimating-equation system (one observation = one id):
  1. Cause-1 hazard (beta1), summed over each id's at-risk PP rows:
       sum_k atrisk_ik * X_ik' * (1{event_ik == 1} - expit(X_ik beta1))
  2. Cause-2 hazard (beta2): same with 1{event_ik == 2}.
     X_ik = [1, factor(t) dummies, A_i, L_i] (the survatr ~ factor(t) + A + L
     design); the shared at-risk set drops rows at/after the first event of ANY
     cause, i.e. competing causes are censored at their event time.
  3. CIF functional (theta_{j,a,t}), per id:
       sum_{k<=t} S^a_i(k-1) h^(j),a_ik - theta_{j,a,t}
  4. All-cause survival functional (theta_{S,a,t}), per id:
       prod_{k<=t} (1 - sum_j h^(j),a_ik) - theta_{S,a,t}

The block-diagonal bread across the two hazard scores plus the cross-time delta
through the CIF / survival functionals is exactly survatr's stacked sandwich, so
the delicatessen variance is an independent check of survatr's analytic IF (in
particular the all-cause (1 - H) sensitivity denominator).

Usage (delicatessen lives in causatr's data-raw venv):
    ../causatr/data-raw/zepid_venv/bin/python delicatessen_competing_risks.py

Reads  tests/testthat/fixtures/python/cr_survival_data.csv  (shared with R)
Writes tests/testthat/fixtures/python/cr_survival_delicatessen.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator

TIMES = [2, 4]  # CIF horizons to pin
CAUSES = [1, 2]


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


def cr_reference(df):
    ids = np.sort(df["id"].unique())
    n = len(ids)
    K = int(df["t"].max())
    p = K + 2  # intercept + (K-1) time dummies + A + L
    periods = np.arange(1, K + 1)

    L = np.array([df.loc[df["id"] == i, "L"].iloc[0] for i in ids])
    A = np.array([df.loc[df["id"] == i, "A"].iloc[0] for i in ids], dtype=float)
    # Per-cause event indicators and the shared at-risk mask (n x K).
    D1 = np.zeros((n, K))
    D2 = np.zeros((n, K))
    atrisk = np.zeros((n, K))
    for r, i in enumerate(ids):
        sub = df[df["id"] == i].sort_values("t")
        ev = sub["event"].values
        D1[r, :] = (ev == 1).astype(float)
        D2[r, :] = (ev == 2).astype(float)
        atrisk[r, :] = sub["at_risk"].values

    Xobs = np.stack([build_design(periods, K, A[r], L[r]) for r in range(n)])
    # Counterfactual designs per id under A = a (keep each id's L).
    X_cf = {
        a: np.stack([build_design(periods, K, float(a), L[r]) for r in range(n)])
        for a in (0, 1)
    }
    Dlist = {1: D1, 2: D2}

    # theta = [beta1(p), beta2(p), CIF(theta_{j,a,t}), S(theta_{a,t})]
    i_b = {1: 0, 2: p}
    n_cif = len(CAUSES) * 2 * len(TIMES)  # cause x arm x time
    i_cif = 2 * p
    i_surv = i_cif + n_cif

    # Index maps for readout.
    cif_idx = {}
    row = 0
    for j in CAUSES:
        for a in (1, 0):
            for t in TIMES:
                cif_idx[(j, a, t)] = i_cif + row
                row += 1
    surv_idx = {}
    row = 0
    for a in (1, 0):
        for t in TIMES:
            surv_idx[(a, t)] = i_surv + row
            row += 1

    def psi(theta):
        b1 = theta[i_b[1]:i_b[1] + p]
        b2 = theta[i_b[2]:i_b[2] + p]
        betas = {1: b1, 2: b2}

        # --- cause-specific hazard scores, summed over at-risk PP rows ---
        ee_b = {}
        for j in CAUSES:
            h = expit(np.einsum("nkp,p->nk", Xobs, betas[j]))
            resid = (Dlist[j] - h) * atrisk
            ee_b[j] = np.einsum("nkp,nk->pn", Xobs, resid)

        # --- per-id counterfactual hazards, all-cause survival, CIF ---
        ee_cif = np.zeros((n_cif, n))
        ee_surv = np.zeros((len(CAUSES) * 0 + 2 * len(TIMES), n))
        for a in (1, 0):
            hj = {j: expit(np.einsum("nkp,p->nk", X_cf[a], betas[j]))
                  for j in CAUSES}
            H = hj[1] + hj[2]  # n x K all-cause hazard
            surv = np.cumprod(1 - H, axis=1)  # S^a_i(k), k=1..K
            # S^a_i(k-1): prepend 1, drop last column.
            surv_lag = np.concatenate([np.ones((n, 1)), surv[:, :-1]], axis=1)
            for j in CAUSES:
                cif_inc = surv_lag * hj[j]  # n x K
                cif_cum = np.cumsum(cif_inc, axis=1)
                for t in TIMES:
                    r = cif_idx[(j, a, t)] - i_cif
                    ee_cif[r, :] = cif_cum[:, t - 1] - theta[cif_idx[(j, a, t)]]
            for t in TIMES:
                r = surv_idx[(a, t)] - i_surv
                ee_surv[r, :] = surv[:, t - 1] - theta[surv_idx[(a, t)]]

        return np.vstack([ee_b[1], ee_b[2], ee_cif, ee_surv])

    # --- initial values: per-cause logistic fits on at-risk rows ---
    def nll(b, Xd, yd):
        pr = np.clip(expit(Xd @ b), 1e-10, 1 - 1e-10)
        return -np.sum(yd * np.log(pr) + (1 - yd) * np.log(1 - pr))

    binit = {}
    Xflat = Xobs.reshape(n * K, p)
    mask = atrisk.reshape(n * K) > 0
    for j in CAUSES:
        yj = Dlist[j].reshape(n * K)
        binit[j] = minimize(
            lambda b: nll(b, Xflat[mask], yj[mask]), np.zeros(p), method="BFGS"
        ).x

    cif_init = np.zeros(n_cif)
    surv_init = np.zeros(2 * len(TIMES))
    for a in (1, 0):
        hj = {j: expit(np.einsum("nkp,p->nk", X_cf[a], binit[j])) for j in CAUSES}
        H = hj[1] + hj[2]
        surv = np.cumprod(1 - H, axis=1)
        surv_lag = np.concatenate([np.ones((n, 1)), surv[:, :-1]], axis=1)
        for j in CAUSES:
            cif_cum = np.cumsum(surv_lag * hj[j], axis=1)
            for t in TIMES:
                cif_init[cif_idx[(j, a, t)] - i_cif] = cif_cum[:, t - 1].mean()
        for t in TIMES:
            surv_init[surv_idx[(a, t)] - i_surv] = surv[:, t - 1].mean()
    init = np.concatenate([binit[1], binit[2], cif_init, surv_init])

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm")
    theta = mest.theta
    vcov = mest.variance

    out = []
    for j in CAUSES:
        for a in (1, 0):
            for t in TIMES:
                k = cif_idx[(j, a, t)]
                out.append(dict(quantity=f"F{j}_a{a}", time=t,
                                estimate=theta[k], se=np.sqrt(vcov[k, k])))
    # All-cause survival.
    for a in (1, 0):
        for t in TIMES:
            k = surv_idx[(a, t)]
            out.append(dict(quantity=f"S_a{a}", time=t,
                            estimate=theta[k], se=np.sqrt(vcov[k, k])))
    # cif_difference for each cause: survatr "a0 vs a1" = F_a0 - F_a1.
    for j in CAUSES:
        for t in TIMES:
            j0, j1 = cif_idx[(j, 0, t)], cif_idx[(j, 1, t)]
            est = theta[j0] - theta[j1]
            var = vcov[j0, j0] + vcov[j1, j1] - 2 * vcov[j0, j1]
            out.append(dict(quantity=f"RD{j}", time=t,
                            estimate=est, se=np.sqrt(var)))
    return pd.DataFrame(out)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fx = os.path.normpath(
        os.path.join(here, "..", "tests", "testthat", "fixtures", "python")
    )
    df = pd.read_csv(os.path.join(fx, "cr_survival_data.csv"))
    print(f"Loaded {len(df)} rows, {df['id'].nunique()} ids, K={df['t'].max()}")
    res = cr_reference(df)
    res["estimate"] = res["estimate"].round(8)
    res["se"] = res["se"].round(8)
    print(res.to_string(index=False))
    res.to_csv(os.path.join(fx, "cr_survival_delicatessen.csv"), index=False)
    print("wrote cr_survival_delicatessen.csv")
