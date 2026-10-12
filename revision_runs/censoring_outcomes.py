"""Censoring aware outcomes for CausalPFN (run D), computed on the training partition only.

Usage: python censoring_outcomes.py IN_CSV OUT_CSV METHOD HORIZON
  IN_CSV  columns patientID, time, event, treatment (training partition)
  OUT_CSV columns patientID, Y
  METHOD  ipcw         Y = Delta_h * min(U, h) / G(min(U, h)-), where Delta_h = 1 if min(U, h) is observed
                       (event before h, or follow up reaching h) and G is the Kaplan Meier survival of the
                       censoring time fitted separately in each treatment arm. G is clipped at 0.05, as in
                       policy_value_analysis/policy_value.py.
          pseudo_rmst  jackknife pseudo observations of RMST(h): Y_i = n * theta - (n - 1) * theta_(-i),
                       theta = integral_0^h of the Kaplan Meier survival of the whole training partition.
Writes a JSON summary next to OUT_CSV.
"""
import json
import sys

import numpy as np
import pandas as pd

G_CLIP = 0.05


def km_curve(time, event):
    """Kaplan Meier S(t) as step function: returns (event times, S just after each)."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    ut = np.unique(t[e == 1])
    at_risk = np.array([(t >= u).sum() for u in ut])
    deaths = np.array([((t == u) & (e == 1)).sum() for u in ut])
    return ut, np.cumprod(1 - deaths / at_risk)


def km_left_limit(ut, S, x):
    """S(x-) for a right continuous step function with jumps at ut."""
    idx = np.searchsorted(ut, x, side="left") - 1
    return np.where(idx < 0, 1.0, S[np.clip(idx, 0, None)]) if len(ut) else np.ones_like(x, dtype=float)


def rmst(time, event, h):
    ut, S = km_curve(time, event)
    keep = ut <= h
    ut, S = ut[keep], S[keep]
    knots = np.concatenate([[0.0], ut, [h]])
    levels = np.concatenate([[1.0], S])
    return float(np.sum(levels * np.diff(knots)))


def ipcw(d, h):
    U = np.minimum(d.time.values, h)
    observed = ((d.event.values == 1) & (d.time.values <= h)) | (d.time.values >= h)
    G = np.ones(len(d))
    for arm in (0, 1):
        g = d.treatment.values == arm
        ut, S = km_curve(d.time.values[g], 1 - d.event.values[g])
        G[g] = km_left_limit(ut, S, U[g])
    n_clipped = int(((G < G_CLIP) & observed).sum())
    G = np.clip(G, G_CLIP, 1.0)
    Y = observed * U / G
    info = {"n_observed_restricted": int(observed.sum()), "n_zero_weight": int((~observed).sum()),
            "min_G_observed": float(G[observed].min()), "n_G_clipped": n_clipped}
    return Y, info


def pseudo_rmst(d, h):
    t, e = d.time.values, d.event.values
    n = len(t)
    theta = rmst(t, e, h)
    loo = np.empty(n)
    mask = np.ones(n, bool)
    for i in range(n):
        mask[i] = False
        loo[i] = rmst(t[mask], e[mask], h)
        mask[i] = True
    Y = n * theta - (n - 1) * loo
    return Y, {"theta_rmst_h": theta, "mean_pseudo": float(Y.mean())}


def main():
    in_csv, out_csv, method, h = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    d = pd.read_csv(in_csv)
    if method == "ipcw":
        Y, info = ipcw(d, h)
    elif method == "pseudo_rmst":
        Y, info = pseudo_rmst(d, h)
    else:
        sys.exit(f"unknown method {method}")
    pd.DataFrame({"patientID": d.patientID, "Y": Y}).to_csv(out_csv, index=False)
    info.update({"method": method, "horizon": h, "n": len(d), "Y_min": float(Y.min()), "Y_median": float(np.median(Y)),
                 "Y_mean": float(Y.mean()), "Y_max": float(Y.max())})
    with open(out_csv.replace(".csv", ".json"), "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info))


if __name__ == "__main__":
    main()
