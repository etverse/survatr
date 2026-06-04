"""
Generate Track B (ICE-hazard survival) reference values via delicatessen
stacked M-estimation, as an INDEPENDENT cross-check of survatr's analytic
survival influence-function sandwich.

This implements the SAME estimator survatr's `estimator = "ice"` does: a
backward iterated-conditional-expectation on the discrete-time hazard with a
survival-tail pseudo-outcome, full treatment/confounder history, and a
counterfactual standardisation mean per static strategy. delicatessen forms the
sandwich numerically from the stacked estimating equations, so agreement to
solver tolerance validates survatr's hand-derived chain -- in particular the
(1 - D_k) failure carry-forward factor, which causatr's plain ICE chain lacks.

DGP: `sim_ice_feedback()` (helper-ice-survival-oracle.R), K = 3 periods,
treatment-confounder feedback. Per-period model design (full history, matching
`ice_build_formula`):
    X1 = [1, A1, L1]
    X2 = [1, A2, L2, A1, L1]            # lag1_A, lag1_L
    X3 = [1, A3, L3, A2, L2, A1, L1]    # lag1_*, lag2_*
Counterfactual designs set the CURRENT period's treatment to a*, leaving lag
columns at observed values (the ICE contract).

Estimand at horizon tau: R^a(tau) = E[1 - S^a(tau)] via the backward recursion;
S^a(tau) = 1 - R^a(tau). RD(tau) = S^1(tau) - S^0(tau) = mu0 - mu1.

Usage:
    cd survatr/data-raw
    # data fixture is generated in R (see fixtures/python/README.md)
    ../../causatr/data-raw/zepid_venv/bin/python delicatessen_ice_survival.py
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator


def _nll_logistic(b, x, y, w):
    """Weighted (at-risk-masked) logistic NLL for sane starting values."""
    p = np.clip(expit(x @ b), 1e-10, 1 - 1e-10)
    return -np.sum(w * (y * np.log(p) + (1 - y) * np.log(1 - p)))


def _init_logit(x, y, w):
    p = x.shape[1]
    res = minimize(_nll_logistic, np.zeros(p), args=(x, y, w), method="BFGS")
    return res.x


def ice_survival_horizon(d, tau):
    """
    Stacked-EE ICE survival at horizon `tau` for static a=1 and a=0 jointly.

    Returns dict with S1, S0 (survival) and RD = S1 - S0, each with its
    sandwich SE.

    theta layout:
      [final beta_tau]          (shared; response = observed event at tau)
      [a=1: beta_1 ... beta_{tau-1}]
      [a=0: beta_1 ... beta_{tau-1}]
      [mu1, mu0, rd]
    """
    n = d["id"].nunique()
    # Wide per-period arrays (id-sorted, one row per id-period -> columns).
    w = d.pivot(index="id", columns="t", values=["A", "L", "Y"]).sort_index()
    A = {k: w[("A", k)].to_numpy(float) for k in range(1, 4)}
    L = {k: w[("L", k)].to_numpy(float) for k in range(1, 4)}
    Y = {k: w[("Y", k)].to_numpy(float) for k in range(1, 4)}
    one = np.ones(n)

    # Observed per-period designs (full history).
    def Xobs(k):
        if k == 1:
            return np.column_stack([one, A[1], L[1]])
        if k == 2:
            return np.column_stack([one, A[2], L[2], A[1], L[1]])
        return np.column_stack([one, A[3], L[3], A[2], L[2], A[1], L[1]])

    # Counterfactual designs: current-period A -> a*, lags observed.
    def Xcf(k, a):
        if k == 1:
            return np.column_stack([one, a * one, L[1]])
        if k == 2:
            return np.column_stack([one, a * one, L[2], A[1], L[1]])
        return np.column_stack([one, a * one, L[3], A[2], L[2], A[1], L[1]])

    # At-risk masks: r_k = 1 if at risk entering period k (no prior event).
    r = {
        1: np.ones(n),
        2: (Y[1] == 0).astype(float),
        3: ((Y[1] == 0) & (Y[2] == 0)).astype(float),
    }

    p_fin = Xobs(tau).shape[1]
    # parameter-block sizes for the backward (pseudo) steps 1..tau-1
    pseudo_p = [Xobs(k).shape[1] for k in range(1, tau)]
    n_pseudo = sum(pseudo_p)

    # offsets
    o_fin = 0
    o_a1 = p_fin
    o_a0 = p_fin + n_pseudo
    o_mu1 = p_fin + 2 * n_pseudo
    o_mu0 = o_mu1 + 1
    o_rd = o_mu1 + 2
    n_theta = o_rd + 1

    def split_pseudo(theta, base):
        """Return list of beta blocks for steps 1..tau-1 from a flat slice."""
        out = []
        off = base
        for pk in pseudo_p:
            out.append(theta[off:off + pk])
            off += pk
        return out

    def psi(theta):
        beta_fin = theta[o_fin:o_fin + p_fin]
        b1 = split_pseudo(theta, o_a1)  # a=1 step betas (1..tau-1)
        b0 = split_pseudo(theta, o_a0)  # a=0 step betas
        mu1 = theta[o_mu1]
        mu0 = theta[o_mu0]
        rd = theta[o_rd]

        rows = []

        # -- final step beta_tau: observed event at tau among at-risk-at-tau --
        Xt = Xobs(tau)
        resid_fin = r[tau] * (Y[tau] - expit(Xt @ beta_fin))
        rows.append(Xt.T * resid_fin)

        # -- backward pseudo steps for each strategy --
        for a, bb in ((1, b1), (0, b0)):
            # Walk k = tau-1 down to 1; q_{k+1} uses next step's beta (or the
            # shared final beta when k+1 == tau).
            for k in range(tau - 1, 0, -1):
                if k + 1 == tau:
                    q_next = expit(Xcf(k + 1, a) @ beta_fin)
                else:
                    q_next = expit(Xcf(k + 1, a) @ bb[k])  # bb[k] = beta_{k+1}
                y_tilde = Y[k] + (1 - Y[k]) * q_next
                Xk = Xobs(k)
                resid = r[k] * (y_tilde - expit(Xk @ bb[k - 1]))
                rows.append(Xk.T * resid)

        # -- standardisation means: mu^a = E[ expit(Xcf(1,a) beta^a_1) ] --
        q1_1 = expit(Xcf(1, 1) @ b1[0])
        q1_0 = expit(Xcf(1, 0) @ b0[0])
        rows.append((q1_1 - mu1)[np.newaxis, :])
        rows.append((q1_0 - mu0)[np.newaxis, :])
        # RD = mu0 - mu1 = S1 - S0 (constant EE row, broadcast to n columns)
        rows.append((np.ones(n) * ((mu0 - mu1) - rd))[np.newaxis, :])

        return np.vstack(rows)

    # -- starting values: fit each block sequentially with plug-in pseudo --
    init = np.zeros(n_theta)
    beta_fin0 = _init_logit(Xobs(tau), Y[tau], r[tau])
    init[o_fin:o_fin + p_fin] = beta_fin0

    for a, off in ((1, o_a1), (0, o_a0)):
        betas = {}
        # build from k = tau-1 down to 1; precompute block offsets
        block_off = []
        oo = off
        for pk in pseudo_p:
            block_off.append(oo)
            oo += pk
        for k in range(tau - 1, 0, -1):
            if k + 1 == tau:
                q_next = expit(Xcf(k + 1, a) @ beta_fin0)
            else:
                q_next = expit(Xcf(k + 1, a) @ betas[k + 1])
            y_tilde = Y[k] + (1 - Y[k]) * q_next
            bk = _init_logit(Xobs(k), np.clip(y_tilde, 1e-6, 1 - 1e-6), r[k])
            betas[k] = bk
            init[block_off[k - 1]:block_off[k - 1] + pseudo_p[k - 1]] = bk
        mu_a = np.mean(expit(Xcf(1, a) @ betas[1]))
        if a == 1:
            init[o_mu1] = mu_a
        else:
            init[o_mu0] = mu_a
    init[o_rd] = init[o_mu0] - init[o_mu1]

    mest = MEstimator(psi, init=init)
    mest.estimate(solver="lm", maxiter=5000)
    theta = mest.theta
    vcov = mest.variance

    mu1, mu0 = theta[o_mu1], theta[o_mu0]
    se_mu1 = np.sqrt(vcov[o_mu1, o_mu1])
    se_mu0 = np.sqrt(vcov[o_mu0, o_mu0])
    rd = theta[o_rd]
    se_rd = np.sqrt(vcov[o_rd, o_rd])
    # S = 1 - mu, so SE(S) = SE(mu).
    return {
        "S1": (1 - mu1, se_mu1),
        "S0": (1 - mu0, se_mu0),
        "RD": (rd, se_rd),
    }


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fx = os.path.join(
        here, "..", "tests", "testthat", "fixtures", "python",
        "ice_survival_data.csv"
    )
    d = pd.read_csv(fx)

    out = []
    for tau in (1, 2, 3):
        if tau == 1:
            # Horizon 1: only the final hazard model; mu^a = E[expit(Xcf(1,a))].
            n = d["id"].nunique()
            w = d[d["t"] == 1].sort_values("id")
            A1 = w["A"].to_numpy(float)
            L1 = w["L"].to_numpy(float)
            Y1 = w["Y"].to_numpy(float)
            one = np.ones(n)
            X1 = np.column_stack([one, A1, L1])

            def psi1(theta):
                beta = theta[:3]
                mu1 = theta[3]
                mu0 = theta[4]
                rd = theta[5]
                resid = Y1 - expit(X1 @ beta)
                q1 = expit(np.column_stack([one, one, L1]) @ beta)
                q0 = expit(np.column_stack([one, 0 * one, L1]) @ beta)
                return np.vstack([
                    X1.T * resid,
                    (q1 - mu1)[np.newaxis, :],
                    (q0 - mu0)[np.newaxis, :],
                    (one * ((mu0 - mu1) - rd))[np.newaxis, :],
                ])

            b0 = _init_logit(X1, Y1, one)
            m1 = np.mean(expit(np.column_stack([one, one, L1]) @ b0))
            m0 = np.mean(expit(np.column_stack([one, 0 * one, L1]) @ b0))
            init = np.concatenate([b0, [m1, m0, m0 - m1]])
            mest = MEstimator(psi1, init=init)
            mest.estimate(solver="lm")
            th, vc = mest.theta, mest.variance
            res = {
                "S1": (1 - th[3], np.sqrt(vc[3, 3])),
                "S0": (1 - th[4], np.sqrt(vc[4, 4])),
                "RD": (th[5], np.sqrt(vc[5, 5])),
            }
        else:
            res = ice_survival_horizon(d, tau)
        for q in ("S1", "S0", "RD"):
            est, se = res[q]
            out.append({"quantity": q, "time": tau, "estimate": est, "se": se})
            print(f"  t={tau} {q}: est={est:.5f} se={se:.5f}")

    ref = pd.DataFrame(out)
    out_path = os.path.join(
        here, "..", "tests", "testthat", "fixtures", "python",
        "ice_survival_delicatessen.csv"
    )
    ref.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}")
