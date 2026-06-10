"""
Reference values for the IPCW IPW weighted-hazard MSM via delicatessen stacked
M-estimation (Zivich 2022). Independent analytic-sandwich oracle for survatr's
`estimator = "ipw"` + `ipcw = ~L` + `ci_method = "sandwich"` path (chunk 11).

The censoring is informative: C_k ~ Bernoulli(expit(cens0 + delta_cens * L)),
where delta_cens = 0.8. Ignoring censoring over-estimates survival; IPCW
corrects by reweighting each at-risk row by W^C_{i,k}.

Stacked estimating-equation system (one observation = one id):
  1. Treatment propensity (alpha, p_alpha):
       Z_i' (A_i - expit(Z_i alpha)),  Z = [1, L_i]
  2. Censoring hazard denominator (gamma, p_gamma):
       sum_{cens_at_risk_ik} X_cens_{ik}' (C_{ik} - expit(X_cens_{ik} gamma))
       X_cens = [1, I(t=2..5), A, L]  (factor(t) + A + L)
  3. Weighted hazard MSM (beta, p_beta), at at-risk rows:
       sum_{at_risk_ik} X_{ik}' * w_T_i(alpha) * W^C_{i,k}(gamma) *
           (Y_{ik} - expit(X_{ik} beta))
       X = [1, I(t=2..5), A]  (factor(t) + A, the survatr MSM design)
  4. Survival functional (theta_{a,t}), per id:
       prod_{k<=t}(1 - expit(X^{A=a}_{ik} beta)) - theta_{a,t}

The stabilizing numerator for both W_T and W^C is treated as FIXED at the
estimated MLE (conservative, consistent with survatr's omission of the
numerator estimation from the IF).

Running IPCW product:
  W^C_{i,k} = prod_{m=1}^{k} (1 - g_num_{im}) / (1 - g_{im}(gamma))
Computed over ALL person-period rows (the at-risk restriction applies only to
which rows enter the MSM score, not which periods contribute to the product).

Usage (delicatessen lives in causatr's data-raw venv):
    ../causatr/data-raw/zepid_venv/bin/python delicatessen_ipcw_survival.py

Reads  tests/testthat/fixtures/python/ipcw_survival_data.csv  (shared with R)
Writes tests/testthat/fixtures/python/ipcw_survival_delicatessen.csv
"""

import os

import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.optimize import minimize
from delicatessen import MEstimator

TIMES = [2, 3, 4, 5]  # survival horizons to pin


def build_factor_t_dummies(periods, K):
    """Columns I(t=2), ..., I(t=K) for each row in `periods`."""
    n = len(periods)
    dummies = np.zeros((n, K - 1))
    for j in range(2, K + 1):
        dummies[:, j - 2] = (periods == j).astype(float)
    return dummies


def ipcw_survival_reference(df):
    ids = np.sort(df["id"].unique())
    n = len(ids)
    K = int(df["t"].max())
    periods = np.arange(1, K + 1)

    # ── Reshape to (n, K) arrays ──────────────────────────────────────────────
    id_to_row = {i: r for r, i in enumerate(ids)}
    A_arr = np.zeros(n)
    L_arr = np.zeros(n)
    Y_arr = np.zeros((n, K))
    C_arr = np.zeros((n, K))
    at_risk_arr = np.zeros((n, K))
    cens_at_risk_arr = np.zeros((n, K))

    for i in ids:
        r = id_to_row[i]
        sub = df[df["id"] == i].sort_values("t")
        A_arr[r] = sub["A"].iloc[0]
        L_arr[r] = sub["L"].iloc[0]
        Y_arr[r, :] = sub["Y"].values
        C_arr[r, :] = sub["C"].values
        at_risk_arr[r, :] = sub["at_risk"].values
        cens_at_risk_arr[r, :] = sub["cens_at_risk"].values

    # ── Design matrices (n, K, p) ─────────────────────────────────────────────
    # Time dummies shape: (K, K-1) then broadcast to (n, K, K-1)
    t_dummies_K = build_factor_t_dummies(periods, K)  # (K, K-1)
    t_dummies_nK = np.broadcast_to(t_dummies_K[np.newaxis, :, :], (n, K, K - 1))

    # Cens denominator: [1, factor(t), A, L] -- p_gamma = K+2 = 7
    p_gamma = K + 2
    X_cens = np.zeros((n, K, p_gamma))
    X_cens[:, :, 0] = 1.0                       # intercept
    X_cens[:, :, 1:K] = t_dummies_nK            # time dummies (4 columns)
    X_cens[:, :, K] = A_arr[:, np.newaxis]       # treatment
    X_cens[:, :, K + 1] = L_arr[:, np.newaxis]  # confounder

    # MSM / numerator: [1, factor(t), A] -- p_beta = p_num = K+1 = 6
    p_beta = K + 1
    X_msm = np.zeros((n, K, p_beta))
    X_msm[:, :, 0] = 1.0
    X_msm[:, :, 1:K] = t_dummies_nK
    X_msm[:, :, K] = A_arr[:, np.newaxis]

    # Counterfactual designs (K, p_beta): A forced to a constant
    X_cf = {}
    for a in (0, 1):
        X_a = np.zeros((K, p_beta))
        X_a[:, 0] = 1.0
        X_a[:, 1:K] = t_dummies_K
        X_a[:, K] = float(a)
        X_cf[a] = X_a

    # Treatment propensity design: [1, L] -- p_alpha = 2
    p_alpha = 2
    Z = np.column_stack([np.ones(n), L_arr])  # (n, p_alpha)

    # ── Numerator models (fixed at MLE; not in stacked EE) ───────────────────
    # Treatment numerator: fixed marginal P(A=1)
    pm_T = A_arr.mean()
    f_num_T = np.where(A_arr == 1, pm_T, 1 - pm_T)  # n

    # Censoring numerator: X_num = [1, factor(t), A] = same layout as X_msm
    # Pre-fit via logistic regression of C on X_num at cens_at_risk rows
    X_num_flat = X_msm.reshape(n * K, p_beta)
    C_flat = C_arr.reshape(n * K)
    mask_cens_fit = cens_at_risk_arr.reshape(n * K).astype(bool)
    from scipy.optimize import minimize as _minimize

    def _nll_num(b):
        p = np.clip(expit(X_num_flat[mask_cens_fit] @ b), 1e-10, 1 - 1e-10)
        c = C_flat[mask_cens_fit]
        return -np.sum(c * np.log(p) + (1 - c) * np.log(1 - p))

    gamma_num_hat = _minimize(_nll_num, np.zeros(p_beta), method="BFGS").x
    # g_num for ALL (n, K) rows (used as fixed factor in W^C)
    eta_num_all = np.einsum("nkp,p->nk", X_msm, gamma_num_hat)  # (n, K)
    g_num_all = expit(eta_num_all)  # (n, K)
    # Fixed numerator factor: (1 - g_num)
    numer_factor = 1.0 - g_num_all  # (n, K)

    # ── Stacked EE ────────────────────────────────────────────────────────────
    n_surv = 2 * len(TIMES)
    i_gamma = p_alpha
    i_beta = p_alpha + p_gamma
    i_surv = p_alpha + p_gamma + p_beta

    def psi(theta):
        alpha = theta[:p_alpha]
        gamma = theta[i_gamma:i_gamma + p_gamma]
        beta = theta[i_beta:i_beta + p_beta]
        surv = theta[i_surv:]

        # ── Block 1: treatment propensity ─────────────────────────────────────
        pi_T = expit(Z @ alpha)                       # n
        ee_alpha = Z.T * (A_arr - pi_T)               # (p_alpha, n)

        # Stabilized treatment weight W_T = f_num_T / P(A | L)
        f_den_T = np.where(A_arr == 1, pi_T, 1 - pi_T)
        w_T = f_num_T / f_den_T                       # n

        # ── Block 2: censoring denominator model ──────────────────────────────
        # g_den = expit(X_cens gamma) for all (n, K)
        eta_cens = np.einsum("nkp,p->nk", X_cens, gamma)  # (n, K)
        g_den = expit(eta_cens)                            # (n, K)
        # score = (C - g_den) at cens_at_risk rows
        score_cens = cens_at_risk_arr * (C_arr - g_den)   # (n, K)
        # per-id sum: X_cens' * score_cens summed over k
        # einsum: nkp,nk -> pn
        ee_gamma = np.einsum("nkp,nk->pn", X_cens, score_cens)  # (p_gamma, n)

        # ── IPCW running product W^C (n, K) ───────────────────────────────────
        # per-row factor = numer_factor / (1 - g_den)  [fixed num, gamma-varying den]
        per_row_factor = numer_factor / (1.0 - g_den)  # (n, K)
        # Running cumprod over time axis (k), broadcasting the product left-to-right
        W_C = np.cumprod(per_row_factor, axis=1)        # (n, K)

        # ── Block 3: weighted hazard MSM ──────────────────────────────────────
        h_obs = expit(np.einsum("nkp,p->nk", X_msm, beta))   # (n, K)
        W_comb = w_T[:, np.newaxis] * W_C                     # (n, K)
        # score at at-risk rows (C_k=0, still in risk set)
        score_msm = at_risk_arr * W_comb * (Y_arr - h_obs)    # (n, K)
        ee_beta = np.einsum("nkp,nk->pn", X_msm, score_msm)  # (p_beta, n)

        # ── Block 4: survival functional ──────────────────────────────────────
        ee_surv = np.zeros((n_surv, n))
        row = 0
        for a in (1, 0):
            # Counterfactual hazard: (K,) -- same for all ids (marginal MSM)
            h_a = expit(X_cf[a] @ beta)       # K
            S_cum = np.cumprod(1.0 - h_a)     # K: S^a(t) = prod_{k<=t}(1-h_a_k)
            for t in TIMES:
                # Per-id EE: S^a_i(t) - theta_{a,t} (simplified: id-invariant)
                ee_surv[row, :] = S_cum[t - 1] - surv[row]
                row += 1

        return np.vstack([ee_alpha, ee_gamma, ee_beta, ee_surv])

    # ── Initial values ────────────────────────────────────────────────────────
    def _nll(b, Xd, yd, wd=None):
        p = np.clip(expit(Xd @ b), 1e-10, 1 - 1e-10)
        w = wd if wd is not None else np.ones(len(yd))
        return -np.sum(w * (yd * np.log(p) + (1 - yd) * np.log(1 - p)))

    # alpha
    a_init = minimize(lambda b: _nll(b, Z, A_arr), np.zeros(p_alpha), method="BFGS").x
    pi0 = expit(Z @ a_init)
    w_T0 = f_num_T / np.where(A_arr == 1, pi0, 1 - pi0)

    # gamma (censoring denominator)
    X_cens_flat = X_cens.reshape(n * K, p_gamma)
    mask_cf = cens_at_risk_arr.reshape(n * K).astype(bool)
    g_init = minimize(
        lambda b: _nll(b, X_cens_flat[mask_cf], C_flat[mask_cf]),
        np.zeros(p_gamma),
        method="BFGS"
    ).x

    # W^C initial
    g_den0 = expit(np.einsum("nkp,p->nk", X_cens, g_init))
    W_C0 = np.cumprod(numer_factor / (1.0 - g_den0), axis=1)

    # beta
    X_msm_flat = X_msm.reshape(n * K, p_beta)
    Y_flat = Y_arr.reshape(n * K)
    W_comb0_flat = (w_T0[:, np.newaxis] * W_C0).reshape(n * K)
    mask_msm = at_risk_arr.reshape(n * K).astype(bool)
    b_init = minimize(
        lambda b: _nll(b, X_msm_flat[mask_msm], Y_flat[mask_msm],
                       W_comb0_flat[mask_msm]),
        np.zeros(p_beta),
        method="BFGS"
    ).x

    # theta
    surv_init = []
    for a in (1, 0):
        Sc = np.cumprod(1.0 - expit(X_cf[a] @ b_init))
        surv_init += [Sc[t - 1] for t in TIMES]
    init = np.concatenate([a_init, g_init, b_init, surv_init])

    # ── Fit ───────────────────────────────────────────────────────────────────
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
            out.append(dict(quantity=f"S{a}", time=t,
                            estimate=theta[j], se=np.sqrt(vcov[j, j])))
    for t in TIMES:
        j1, j0 = surv_idx[(1, t)], surv_idx[(0, t)]
        est = theta[j1] - theta[j0]
        var_rd = vcov[j1, j1] + vcov[j0, j0] - 2 * vcov[j1, j0]
        out.append(dict(quantity="RD", time=t, estimate=est, se=np.sqrt(var_rd)))
    return pd.DataFrame(out)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    fx = os.path.normpath(
        os.path.join(here, "..", "tests", "testthat", "fixtures", "python")
    )
    df = pd.read_csv(os.path.join(fx, "ipcw_survival_data.csv"))
    print(f"Loaded {len(df)} rows, {df['id'].nunique()} ids, K={df['t'].max()}")
    res = ipcw_survival_reference(df)
    res["estimate"] = res["estimate"].round(8)
    res["se"] = res["se"].round(8)
    print(res.to_string(index=False))
    res.to_csv(os.path.join(fx, "ipcw_survival_delicatessen.csv"), index=False)
    print("wrote ipcw_survival_delicatessen.csv")
