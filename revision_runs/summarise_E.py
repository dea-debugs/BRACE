"""Aggregate run E (22 inputs x 7 seeds) into E_seed_summary.csv and E_pam50_paired.csv.

E_seed_summary.csv: per input, mean and 2.5/97.5 percentiles across the seven seeds (numpy linear
interpolation, so with 7 values they sit just inside min and max) of the pooled, no chemo arm and chemo arm
RMST differences at 150 months and of the number recommended chemo; ranked by mean pooled RMST.
Includes a reference row for the policy "no chemo for anyone".
E_pam50_paired.csv: per input, the per seed (same split) difference PAM50 clean minus input in pooled RMST.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from run_EG import INPUTS, SEEDS

REV = Path("/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs")
PAPER_TABLE3 = {"pearson100": 16.74, "gs1_driver_core": 17.92, "gs2_driver_ext_wgs": 13.51, "gs3_driver_found": 15.43,
                "gs4A_tcga_cohorts": 12.79, "gs4B_tcga_subgroups": 17.12, "gs4C_tcga_all": 15.92,
                "gs5_interval_screen": 17.23, "gs6_prolif_core_cc": 14.52, "gs7_pam50_clean": 18.94,
                "gs8_mammaprint": 16.89, "gs9_blueprint": 17.91, "gs10_oncotype": 8.96, "gs11_endopredict": 17.12,
                "gs12_index_bci": 16.47, "ct_I_k": 17.12, "ct_I_kin": 18.77, "ct_I_kout": 16.55, "ct_II_k": 16.25,
                "ct_II_kin": 16.73, "ct_II_kout": 17.67}

rows, long = [], []
for key, desc, _ in INPUTS:
    for s in SEEDS:
        label = f"E_{key}_s{s}"
        ev = json.loads((REV / f"{label}_eval.json").read_text())
        cfg = json.loads((REV / f"{label}.json").read_text())
        long.append({"key": key, "input": desc, "seed": s, "n_features": cfg["n_features"],
                     "pooled": ev["pooled"]["diff"], "pooled_p": ev["pooled"]["p"],
                     "no_chemo_arm": ev["no_chemo_arm"]["diff"], "chemo_arm": ev["chemo_arm"]["diff"],
                     "n_rec_chemo": ev["n_rec_chemo"],
                     "n_test_chemo": ev["chemo_arm"]["n"], "n_test_no_chemo": ev["no_chemo_arm"]["n"]})
L = pd.DataFrame(long)
L.to_csv(REV / "E_seed_long.csv", index=False)


def q(x, p):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return float(np.percentile(x, p)) if len(x) else np.nan


for key, g in L.groupby("key", sort=False):
    r = {"input": g.input.iloc[0], "key": key, "n_features": int(g.n_features.iloc[0]), "n_seeds": len(g)}
    for m in ("pooled", "no_chemo_arm", "chemo_arm", "n_rec_chemo"):
        r[f"{m}_mean"] = float(g[m].mean())
        r[f"{m}_p2.5"] = q(g[m], 2.5)
        r[f"{m}_p97.5"] = q(g[m], 97.5)
    r["n_seeds_pooled_p_lt_0.05"] = int((g.pooled_p < 0.05).sum())
    r["n_seeds_zero_rec_chemo"] = int((g.n_rec_chemo == 0).sum())  # degenerate: identical to the reference row
    r["seed123_pooled"] = float(g.loc[g.seed == 123, "pooled"].iloc[0])
    r["paper_table3_pooled"] = PAPER_TABLE3.get(key, np.nan)
    rows.append(r)
# Reference policy "no chemo for anyone": Followed = untreated, Not Followed = treated. Taken per seed from
# runs with zero recommended chemo (identical for every such run); within arm differences are undefined.
deg = L[L.n_rec_chemo == 0].groupby("seed").pooled
assert (deg.max() - deg.min()).abs().max() < 1e-9 and set(deg.groups) == set(SEEDS), "need a degenerate run per seed"
ref = deg.first().loc[SEEDS]
rows.append({"input": "REFERENCE: no chemo for anyone", "key": "reference_treat_none", "n_features": 0,
             "n_seeds": len(ref), "pooled_mean": float(ref.mean()), "pooled_p2.5": q(ref, 2.5),
             "pooled_p97.5": q(ref, 97.5), "n_rec_chemo_mean": 0.0, "n_rec_chemo_p2.5": 0.0, "n_rec_chemo_p97.5": 0.0,
             "seed123_pooled": float(ref.loc[123])})
S = pd.DataFrame(rows).sort_values("pooled_mean", ascending=False).reset_index(drop=True)
S.insert(0, "rank", np.arange(1, len(S) + 1))
S.to_csv(REV / "E_seed_summary.csv", index=False)

# paired comparison: PAM50 clean minus each input, same split per seed
P = L.pivot(index="seed", columns="key", values="pooled")
pam = P["gs7_pam50_clean"]
paired = []
for key, desc, _ in INPUTS:
    if key == "gs7_pam50_clean":
        continue
    d = pam - P[key]
    paired.append({"input": desc, "key": key, "mean_diff_pam50_minus_input": float(d.mean()),
                   "min_diff": float(d.min()), "max_diff": float(d.max()), "n_seeds_pam50_higher": int((d > 0).sum())})
pd.DataFrame(paired).sort_values("mean_diff_pam50_minus_input").to_csv(REV / "E_pam50_paired.csv", index=False)

# seed 150 arm counts (checklist 33)

ev150 = json.loads((REV / "E_gs7_pam50_clean_s150_eval.json").read_text())
print("seed 150 test arms (chemo / no chemo):", ev150["chemo_arm"]["n"], "/", ev150["no_chemo_arm"]["n"])
print(S[["rank", "input", "n_features", "pooled_mean", "pooled_p2.5", "pooled_p97.5", "no_chemo_arm_mean",
         "chemo_arm_mean", "n_rec_chemo_mean", "n_seeds_zero_rec_chemo", "seed123_pooled",
         "paper_table3_pooled"]].round(2).to_string(index=False))
