"""Re-rank run E inputs on metrics with a non trivial baseline -> E_policy_ranking.csv.

Per input, across the seven seeds (each seed = its own split): mean, min and max of
  overlap V(BRACE) - V(observed)   (DR, overlap population 0.1 <= e(x) <= 0.9, nuisance clinical+pam50 on train)
  overlap V(BRACE) - V(treat none) (DR)
  no chemo arm RMST difference at 150 months (Followed - Not Followed; undefined when nobody is Not Followed)
plus the share of seeds where the policy is degenerate (chemo for nobody or for everybody) and the number of
seeds whose bootstrap 95% CI for V(BRACE) - V(observed) lies entirely above zero.
Reference rows use the same seven overlap populations (identical across inputs within a seed):
  no chemo for anyone: V(treat none) DR; chemo for everyone: V(treat all) DR; observed care: V(observed).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_EG import INPUTS, SEEDS

REV = Path("/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs")
S = pd.read_csv(REV / "summary.csv").set_index("run_label")


def stats(prefix, vals):
    v = np.asarray(vals, float)
    ok = v[~np.isnan(v)]
    return {f"{prefix}_mean": ok.mean() if len(ok) else np.nan, f"{prefix}_min": ok.min() if len(ok) else np.nan,
            f"{prefix}_max": ok.max() if len(ok) else np.nan, f"{prefix}_n_seeds": len(ok)}


rows, refs = [], {"none": [], "all": [], "obs": []}
for key, desc, _ in INPUTS:
    labs = [f"E_{key}_s{s}" for s in SEEDS]
    r = S.loc[labs]
    n_rec = r.n_rec_chemo.to_numpy()
    n_tot = (r.n_rec_chemo + r.n_rec_no_chemo).to_numpy()
    row = {"input": desc, "key": key, "n_features": int(r.n_features.iloc[0])}
    row.update(stats("V_minus_observed", r.pv_overlap_BRACE_minus_observed_DR))
    row["n_seeds_V_minus_observed_CI_above_0"] = int((r.pv_overlap_BRACE_minus_observed_DR_lo > 0).sum())
    row.update(stats("V_minus_treat_none", r.pv_overlap_BRACE_minus_treat_none_DR))
    row["n_seeds_V_minus_treat_none_CI_above_0"] = int((r.pv_overlap_BRACE_minus_treat_none_DR_lo > 0).sum())
    row.update(stats("no_chemo_arm_rmst", r.rmst150_no_chemo_arm))
    row["share_seeds_degenerate"] = float(np.mean((n_rec == 0) | (n_rec == n_tot)))
    row["n_rec_chemo_mean"] = float(n_rec.mean())
    rows.append(row)
    for s, lab in zip(SEEDS, labs):  # reference policies from this run's overlap population
        pv = json.loads((REV / f"{lab}_policy_value.json").read_text())["overlap_policy_values"]
        refs["none"].append((s, pv["treat_none"]["DR"], pv["observed"]["IPW"]))
        refs["all"].append((s, pv["treat_all"]["DR"], pv["observed"]["IPW"]))

# references: one value per seed (check identical across inputs)
ref_rows = []
for name, desc in (("none", "REFERENCE: no chemo for anyone"), ("all", "REFERENCE: chemo for everyone")):
    d = pd.DataFrame(refs[name], columns=["seed", "V", "obs"])
    g = d.groupby("seed")
    assert (g.V.max() - g.V.min()).max() < 1e-6 and (g.obs.max() - g.obs.min()).max() < 1e-6
    v, obs = g.V.first().loc[SEEDS], g.obs.first().loc[SEEDS]
    none = pd.DataFrame(refs["none"], columns=["seed", "V", "obs"]).groupby("seed").V.first().loc[SEEDS]
    row = {"input": desc, "key": f"reference_{name}", "n_features": 0}
    row.update(stats("V_minus_observed", v - obs))
    row.update(stats("V_minus_treat_none", v - none))
    row.update(stats("no_chemo_arm_rmst", [np.nan] * len(SEEDS)))
    row["share_seeds_degenerate"] = 1.0
    ref_rows.append(row)
obs = pd.DataFrame(refs["none"], columns=["seed", "V", "obs"]).groupby("seed")
row = {"input": "REFERENCE: observed care", "key": "reference_observed", "n_features": 0}
row.update(stats("V_minus_observed", np.zeros(len(SEEDS))))
row.update(stats("V_minus_treat_none", (obs.obs.first() - obs.V.first()).loc[SEEDS]))
row.update(stats("no_chemo_arm_rmst", [np.nan] * len(SEEDS)))
row["share_seeds_degenerate"] = np.nan
ref_rows.append(row)

R = pd.DataFrame(rows + ref_rows).sort_values("V_minus_observed_mean", ascending=False).reset_index(drop=True)
R.insert(0, "rank", np.arange(1, len(R) + 1))
R.to_csv(REV / "E_policy_ranking.csv", index=False)

show = R[["rank", "input", "V_minus_observed_mean", "V_minus_observed_min", "V_minus_observed_max",
          "n_seeds_V_minus_observed_CI_above_0", "V_minus_treat_none_mean", "V_minus_treat_none_min",
          "V_minus_treat_none_max", "no_chemo_arm_rmst_mean", "no_chemo_arm_rmst_min", "no_chemo_arm_rmst_max",
          "no_chemo_arm_rmst_n_seeds", "share_seeds_degenerate"]]
print(show.round(2).to_string(index=False))
inp = R[~R.key.str.startswith("reference")]
print("\ninputs with V - observed range entirely above 0 (min across seeds > 0):", int((inp.V_minus_observed_min > 0).sum()))
print("inputs with mean V - observed > 0:", int((inp.V_minus_observed_mean > 0).sum()))
print("runs (of 154) whose bootstrap CI for V - observed is entirely above 0:", int(inp.n_seeds_V_minus_observed_CI_above_0.sum()))
print("runs (of 154) whose bootstrap CI for V - treat none is entirely above 0:", int(inp.n_seeds_V_minus_treat_none_CI_above_0.sum()))
