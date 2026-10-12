"""Run G: recommendation intervals and abstention (checklist 17).

Design (agreed): seed 123 split, test set fixed (594 patients); CausalPFN refitted on seven random 80%
subsamples of the 1,386 training patients (subsample seeds = the seven S2 seeds), PAM50 clean + clinical
(C3 specification). Per patient: min, max, SD of est_effect and number of positive estimates.
Abstention rule: recommend only when all seven estimates agree on the sign.

Writes G_patient_spread.csv and G_confident_consensus.csv (+ .json settings) for evaluate_run.R and
policy_value_for_run.py, plus C3_restricted_to_G_confident.csv (C3's own recommendations on the same patients).
"""
import json

import numpy as np
import pandas as pd

from brace_paths import REV
from run_EG import SEEDS


runs = []
for s in SEEDS:
    r = pd.read_csv(REV / f"G_pam50clean_clinical_sub80_s{s}.csv")
    runs.append(r.set_index("patientID")["est_effect"].rename(f"s{s}"))
E = pd.concat(runs, axis=1)
base = pd.read_csv(REV / "C3_pam50_clinical_seed123.csv").set_index("patientID")
assert len(E) == 594 and set(E.index) == set(base.index)
E = E.loc[base.index]

n_pos = (E > 0).sum(axis=1)
spread = pd.DataFrame({
    "treatment": base.treatment, "time": base.time, "event": base.event,
    "est_effect_C3_full_train": base.est_effect,
    "est_min": E.min(axis=1), "est_max": E.max(axis=1), "est_mean": E.mean(axis=1), "est_sd": E.std(axis=1),
    "n_positive_of_7": n_pos, "all_agree": (n_pos == 0) | (n_pos == 7),
    "consensus_recommendation": np.where(n_pos == 7, 1, np.where(n_pos == 0, 0, np.nan)),
}).join(E)
spread.reset_index().to_csv(REV / "G_patient_spread.csv", index=False)

conf = spread[spread.all_agree].reset_index()
label = "G_confident_consensus"
out = pd.DataFrame({"patientID": conf.patientID, "treatment": conf.treatment, "time": conf.time, "event": conf.event,
                    "est_effect": conf.est_mean, "recommendation": conf.consensus_recommendation.astype(int),
                    "run_label": label})
out.to_csv(REV / f"{label}.csv", index=False)
(REV / f"{label}.json").write_text(json.dumps({
    "run_label": label, "seed": 123, "n_features": 65, "outcome": "rst_30",
    "rule": "recommend only when all 7 subsample fits (80% of training, seeds " + ", ".join(map(str, SEEDS))
            + ") agree on the sign of est_effect; others abstain", "n_confident": len(out), "n_test": 594}, indent=2))

# C3's own recommendations on the confident patients, for a like for like comparison
c3r = base.loc[conf.patientID].reset_index()
c3r = c3r.assign(run_label="C3_restricted_to_G_confident")[["patientID", "treatment", "time", "event", "est_effect",
                                                            "recommendation", "run_label"]]
c3r.to_csv(REV / "C3_restricted_to_G_confident.csv", index=False)
(REV / "C3_restricted_to_G_confident.json").write_text(json.dumps({
    "run_label": "C3_restricted_to_G_confident", "seed": 123, "n_features": 65, "outcome": "rst_30"}, indent=2))

c3_rec = base.recommendation
agree_c3 = (spread.loc[spread.all_agree, "consensus_recommendation"] == c3_rec[spread.all_agree]).mean()
info = {
    "n_test": 594, "n_confident": int(spread.all_agree.sum()), "share_confident": float(spread.all_agree.mean()),
    "confident_rec_chemo": int((out.recommendation == 1).sum()), "confident_rec_no_chemo": int((out.recommendation == 0).sum()),
    "abstain": int((~spread.all_agree).sum()),
    "confident_by_observed_arm": {"no_chemo": int(spread.all_agree[spread.treatment == 0].sum()),
                                   "chemo": int(spread.all_agree[spread.treatment == 1].sum())},
    "share_confident_among_C3_rec_chemo": float(spread.all_agree[c3_rec == 1].mean()),
    "share_confident_among_C3_rec_no_chemo": float(spread.all_agree[c3_rec == 0].mean()),
    "consensus_matches_C3_full_train_rec": float(agree_c3),
    "median_range_est_effect": float((spread.est_max - spread.est_min).median()),
    "n_positive_distribution": {int(k): int(v) for k, v in n_pos.value_counts().sort_index().items()},
}
(REV / "G_summary.json").write_text(json.dumps(info, indent=2))
print(json.dumps(info, indent=2))
