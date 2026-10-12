"""Run F: semi synthetic benchmark (checklist 6, Reviewer 2).

Real METABRIC covariates (PAM50 clean 50 genes + 10 clinical with missing indicators, seed 123 split:
1,386 training rows, 594 test rows; imputation fitted on training rows) with simulated treatment and survival.

Simulation (z-scores use training means/SDs):
  propensity   e(x) = expit(a0 + a1 * z), z = standardised (-0.6 age + 0.5 NPI - 0.6 ER + 0.4 ERBB2 - 0.3 FOXA1),
               a1 = 1 (fixed confounding strength), a0 solved so the mean of e(x) over the 1,980 patients equals
               the treated share (50%, 20%, 10%)
  survival     Weibull AFT: log T = mu0(x) + W * tau(x) + 0.8 * log(E), E ~ Exp(1)
               mu0(x) = log(200) - 0.35 NPI + 0.20 ER - 0.15 age - 0.15 MKI67
               tau(x) = -0.40 ESR1 + 0.30 MKI67 - 0.30 age   (log time scale, sign heterogeneous)
  truth        tau_150(x) = RMST_1(x; 150) - RMST_0(x; 150), from the Weibull survival functions
  censoring    C = c_max * Uniform(0, 1), c_max solved to give 0%, 20%, 40%, 60% censored (none for 0%)
  noise        0, 50, 200 independent N(0, 1) features appended
Grid 3 x 4 x 3 = 36 settings, 10 replicates each.

Estimators (fit on the 1,386 training rows, evaluated on the 594 test rows):
  pfn_rst360   CausalPFN (run_causalpfn.py defaults, seed 123) with the paper's outcome min(U, 360)
  pfn_ipcw150  CausalPFN with the run D(i) IPCW pseudo outcome at 150 months (G by arm)
  csf          grf causal_survival_forest, target RMST, horizon 150, 2000 trees
Metrics: PEHE of tau_150, Spearman correlation, sign agreement, true value of the recommended policy
(treat if tau_hat > 0) versus the true optimal policy, treat none and treat all (true RMST at 150).
pfn_rst360 estimates an effect on a different scale (horizon 360), so its PEHE is not comparable; its sign
agreement and policy value are.

Usage: python F_semisynthetic.py [--smoke] [--threads N]
Resumable: a replicate is skipped when its output file exists. Log: revision_runs/logs/F/F_run.log
"""
import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq
from scipy.special import expit

REV = Path("/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs")
OUT = REV / "F_semisynthetic"
LOG_DIR = REV / "logs" / "F"
RSCRIPT_CSF = REV / "R" / "csf_fit.R"
H = 150.0
SIGMA = 0.8
SHARES = [0.5, 0.2, 0.1]
CENSORING = [0.0, 0.2, 0.4, 0.6]
NOISE = [0, 50, 200]
N_REP = 10

sys.path.insert(0, str(REV))
from censoring_outcomes import ipcw  # noqa: E402  (same IPCW outcome as run D)


# ---------------- covariates ----------------
def load_covariates():
    cache = OUT / "covariates.csv"
    if cache.exists():
        return pd.read_csv(cache)
    genes = pd.read_csv(REV / "7_PAM50_clean.csv").Gene.tolist()
    expr = pd.read_csv("/Users/sindypinero/Downloads/4_Deanna/Processed_Datasets/dataset_genes.csv",
                       usecols=["patientID"] + genes)
    clin = pd.read_csv(REV / "clinical_covariates.csv")
    test_ids = set(pd.read_csv(REV / "A_baseline_pam50_seed123.csv").patientID)
    d = expr.merge(clin, on="patientID", validate="one_to_one")
    d.insert(1, "split", np.where(d.patientID.isin(test_ids), "test", "train"))
    assert (d.split == "test").sum() == 594 and (d.split == "train").sum() == 1386
    tr = d.split == "train"
    for c in [c for c in clin.columns if c != "patientID"]:
        if d.loc[tr, c].isna().any():
            d[c + "_miss"] = d[c].isna().astype(float)
        d[c] = d[c].fillna(d.loc[tr, c].median())
    d.to_csv(cache, index=False)
    return d


def zscore(d, col):
    tr = d.split == "train"
    return ((d[col] - d.loc[tr, col].mean()) / d.loc[tr, col].std()).to_numpy()


def weibull_rmst(mu, h=H, n_grid=601):
    t = np.linspace(0, h, n_grid)
    S = np.exp(-(t[None, :] / np.exp(mu)[:, None]) ** (1 / SIGMA))
    return np.trapezoid(S, t, axis=1)


# ---------------- one replicate ----------------
def simulate(d, share, cens, noise, rng):
    age, npi, er = zscore(d, "AGE_AT_DIAGNOSIS"), zscore(d, "NPI"), zscore(d, "ER_IHC")
    esr1, mki67, erbb2, foxa1 = (zscore(d, g) for g in ("ESR1", "MKI67", "ERBB2", "FOXA1"))
    z = -0.6 * age + 0.5 * npi - 0.6 * er + 0.4 * erbb2 - 0.3 * foxa1
    z = (z - z.mean()) / z.std()
    a0 = brentq(lambda a: expit(a + z).mean() - share, -20, 20)
    e = expit(a0 + z)
    W = rng.binomial(1, e)

    mu0 = np.log(200) - 0.35 * npi + 0.20 * er - 0.15 * age - 0.15 * mki67
    tau_log = -0.40 * esr1 + 0.30 * mki67 - 0.30 * age
    R0, R1 = weibull_rmst(mu0), weibull_rmst(mu0 + tau_log)
    logE = np.log(rng.exponential(size=len(d)))
    T = np.exp(mu0 + W * tau_log + SIGMA * logE)

    if cens == 0:
        C = np.full(len(d), np.inf)
    else:
        u = rng.uniform(size=len(d))
        cmax = brentq(lambda c: (c * u < T).mean() - cens, 1e-6, 1e6)
        C = cmax * u
    U = np.minimum(T, C)
    D = (T <= C).astype(int)

    feats = [c for c in d.columns if c not in ("patientID", "split")]
    X = d[feats].to_numpy(float)
    if noise:
        X = np.hstack([X, rng.normal(size=(len(d), noise))])
    return dict(X=X, W=W, U=U, D=D, e=e, R0=R0, R1=R1, tau=R1 - R0)


def fit_pfn(Xtr, Wtr, Ytr, Xte):
    import torch
    from causalpfn import CATEEstimator
    np.random.seed(123); torch.manual_seed(123)
    est = CATEEstimator(device="cpu")
    est.fit(Xtr.astype(np.float32), Wtr.astype(np.float32), Ytr.astype(np.float32))
    return est.estimate_cate(Xte.astype(np.float32))


def fit_csf(Xtr, Wtr, Utr, Dtr, Xte, seed, threads, workdir):
    workdir.mkdir(parents=True, exist_ok=True)
    cols = [f"x{i}" for i in range(Xtr.shape[1])]
    tr = pd.DataFrame(Xtr, columns=cols); tr.insert(0, "W", Wtr); tr.insert(0, "D", Dtr); tr.insert(0, "U", Utr)
    tr.to_csv(workdir / "train.csv", index=False)
    pd.DataFrame(Xte, columns=cols).to_csv(workdir / "test.csv", index=False)
    r = subprocess.run(["Rscript", str(RSCRIPT_CSF), str(workdir / "train.csv"), str(workdir / "test.csv"),
                        str(workdir / "pred.csv"), str(seed), str(H), str(threads)], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return pd.read_csv(workdir / "pred.csv").tau_hat.to_numpy()


def metrics(tau_hat, s, te):
    tau, R0, R1 = s["tau"][te], s["R0"][te], s["R1"][te]
    rec = tau_hat > 0
    v_rec = np.where(rec, R1, R0).mean()
    v_opt = np.maximum(R0, R1).mean()
    return {
        "pehe": float(np.sqrt(np.mean((tau_hat - tau) ** 2))),
        "spearman": float(stats.spearmanr(tau_hat, tau)[0]) if np.std(tau_hat) > 0 else np.nan,
        "sign_agreement": float(np.mean(np.sign(tau_hat) == np.sign(tau))),
        "share_rec_treat": float(rec.mean()),
        "V_recommended": float(v_rec), "V_optimal": float(v_opt), "regret": float(v_opt - v_rec),
        "V_treat_none": float(R0.mean()), "V_treat_all": float(R1.mean()),
    }


def run_replicate(d, k, share, cens, noise, rep, threads, log):
    sid = f"t{int(share * 100):02d}_c{int(cens * 100):02d}_n{noise:03d}"
    out_file = OUT / "reps" / f"{sid}_r{rep:02d}.csv"
    if out_file.exists():
        return None
    rng = np.random.default_rng([2026, k, rep])
    s = simulate(d, share, cens, noise, rng)
    tr, te = (d.split == "train").to_numpy(), (d.split == "test").to_numpy()
    info = {"setting": sid, "treated_share": share, "censoring_target": cens, "noise": noise, "replicate": rep,
            "treated_share_obs": float(s["W"].mean()), "censoring_obs": float(1 - s["D"].mean()),
            "share_true_benefit_test": float((s["tau"][te] > 0).mean()), "n_features": s["X"].shape[1]}
    rows = []
    t0 = time.time()
    y360 = np.minimum(s["U"][tr], 360)
    rows.append({"estimator": "pfn_rst360", **metrics(fit_pfn(s["X"][tr], s["W"][tr], y360, s["X"][te]), s, te),
                 "seconds": time.time() - t0})
    t0 = time.time()
    dtr = pd.DataFrame({"time": s["U"][tr], "event": s["D"][tr], "treatment": s["W"][tr]})
    y_ipcw, _ = ipcw(dtr, H)
    rows.append({"estimator": "pfn_ipcw150", **metrics(fit_pfn(s["X"][tr], s["W"][tr], y_ipcw, s["X"][te]), s, te),
                 "seconds": time.time() - t0})
    t0 = time.time()
    tau_csf = fit_csf(s["X"][tr], s["W"][tr], s["U"][tr], s["D"][tr], s["X"][te], seed=1000 * k + rep,
                      threads=threads, workdir=OUT / "work" / f"{sid}_r{rep:02d}")
    rows.append({"estimator": "csf", **metrics(tau_csf, s, te), "seconds": time.time() - t0})
    res = pd.DataFrame([{**info, **r} for r in rows])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out_file, index=False)
    log(f"{sid} rep {rep}: " + "; ".join(f"{r['estimator']} PEHE {r['pehe']:.1f} sign {r['sign_agreement']:.2f} "
                                         f"regret {r['regret']:.2f} ({r['seconds']:.0f}s)" for r in rows))
    return res


def aggregate():
    files = sorted((OUT / "reps").glob("*.csv"))
    if not files:
        return
    a = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    keys = ["setting", "treated_share", "censoring_target", "noise", "estimator"]
    m = ["pehe", "spearman", "sign_agreement", "share_rec_treat", "V_recommended", "V_optimal", "regret",
         "V_treat_none", "V_treat_all", "censoring_obs", "treated_share_obs", "seconds"]
    g = a.groupby(keys)[m]
    summ = g.mean().add_suffix("_mean").join(g.std().add_suffix("_sd")).join(g.size().rename("n_rep")).reset_index()
    summ.to_csv(OUT / "F_summary.csv", index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true", help="one setting (50%% treated, 40%% censored, 50 noise), one replicate")
    ap.add_argument("--threads", type=int, default=4, help="threads for grf")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True); LOG_DIR.mkdir(parents=True, exist_ok=True)
    logf = open(LOG_DIR / "F_run.log", "a")

    def log(msg):
        line = time.strftime("%Y-%m-%d %H:%M:%S ") + msg
        logf.write(line + "\n"); logf.flush(); print(line, flush=True)

    d = load_covariates()
    grid = list(itertools.product(SHARES, CENSORING, NOISE))
    if args.smoke:
        k = grid.index((0.5, 0.4, 50))
        t0 = time.time()
        res = run_replicate(d, k, 0.5, 0.4, 50, 0, args.threads, log)
        log(f"smoke test: {time.time() - t0:.0f} s" + ("" if res is not None else " (already existed)"))
        if res is not None:
            print(res.T.to_string(header=False))
        return
    log(f"grid start: {len(grid)} settings x {N_REP} replicates")
    for k, (share, cens, noise) in enumerate(grid):
        for rep in range(N_REP):
            try:
                run_replicate(d, k, share, cens, noise, rep, args.threads, log)
            except Exception as ex:  # log and continue; rerun picks up missing replicates
                log(f"FAILED t{share} c{cens} n{noise} rep {rep}: {ex}")
        aggregate()
    aggregate()
    log("grid finished")


if __name__ == "__main__":
    main()
