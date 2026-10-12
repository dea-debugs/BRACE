"""Append one row for a run to revision_runs/summary.csv (created with a header if missing).

Usage: python append_summary.py RUN_LABEL INPUT_DESCRIPTION [REFERENCE_LABEL]
Reads RUN_LABEL.json (settings), RUN_LABEL_eval.json (stratified RMST) and, if present,
RUN_LABEL_policy_value.json (overlap population policy value). REFERENCE_LABEL, if given, is the run
whose recommendations this run is compared with (n_changed_vs_ref). A label already in summary.csv is replaced.
"""
import json
import sys
from pathlib import Path

import pandas as pd

REV = Path("/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs")
SUMMARY = REV / "summary.csv"

OUTCOME_ENCODING = {
    "rst_30": "min(U, 360) (paper encoding)",
    "ipcw150": "IPCW pseudo outcome, h = 150, G by arm (train)",
    "pseudo_rmst150": "jackknife pseudo observation of RMST(150) (train)",
}

label, input_desc = sys.argv[1], sys.argv[2]
ref = sys.argv[3] if len(sys.argv) > 3 else ""

cfg = json.loads((REV / f"{label}.json").read_text())
ev = json.loads((REV / f"{label}_eval.json").read_text())
pv_path = REV / f"{label}_policy_value.json"
pv = json.loads(pv_path.read_text()) if pv_path.exists() else None

row = {
    "run_label": label,
    "input": input_desc,
    "n_features": cfg["n_features"],
    "outcome_encoding": OUTCOME_ENCODING.get(cfg["outcome"], cfg["outcome"]),
    "seed": cfg["seed"],
    "n_rec_chemo": ev["n_rec_chemo"],
    "n_rec_no_chemo": ev["n_rec_no_chemo"],
}
for key, name in (("no_chemo_arm", "no_chemo_arm"), ("chemo_arm", "chemo_arm"), ("pooled", "pooled")):
    g = ev[key]
    row[f"rmst150_{name}"] = g["diff"]
    row[f"rmst150_{name}_lo"] = g["lower"]
    row[f"rmst150_{name}_hi"] = g["upper"]
    row[f"rmst150_{name}_p"] = g["p"]
    row[f"{name}_n_followed_notfollowed"] = f'{g["n_followed"]}/{g["n_not_followed"]}'

if pv:
    pvo = pv["overlap_policy_values"]
    c = pv["bootstrap_CI"]
    row.update({
        "pv_overlap_n": pv["n_overlap"],
        "pv_overlap_V_BRACE_DR": pvo["BRACE"]["DR"],
        "pv_overlap_V_BRACE_DR_lo": c["overlap: V(BRACE) DR"][0],
        "pv_overlap_V_BRACE_DR_hi": c["overlap: V(BRACE) DR"][1],
        "pv_overlap_BRACE_minus_observed_DR": pv["overlap_contrasts"]["V(BRACE) - V(observed) DR"],
        "pv_overlap_BRACE_minus_observed_DR_lo": c["overlap: V(BRACE) - V(observed) DR"][0],
        "pv_overlap_BRACE_minus_observed_DR_hi": c["overlap: V(BRACE) - V(observed) DR"][1],
        "pv_overlap_BRACE_minus_treat_none_DR": pv["overlap_contrasts"]["V(BRACE) - V(treat none) DR"],
        "pv_overlap_BRACE_minus_treat_none_DR_lo": c["overlap: V(BRACE) - V(treat none) DR"][0],
        "pv_overlap_BRACE_minus_treat_none_DR_hi": c["overlap: V(BRACE) - V(treat none) DR"][1],
        "pv_bootstrap_B": pv["bootstrap_B"],
    })

if ref:
    a = pd.read_csv(REV / f"{label}.csv")[["patientID", "recommendation"]]
    b = pd.read_csv(REV / f"{ref}.csv")[["patientID", "recommendation"]]
    m = a.merge(b, on="patientID", suffixes=("", "_ref"), validate="one_to_one")
    assert len(m) == len(a) == len(b)
    row["ref_label"] = ref
    row["n_changed_vs_ref"] = int((m.recommendation != m.recommendation_ref).sum())
    row["changed_0to1_1to0"] = (f'{int(((m.recommendation_ref == 0) & (m.recommendation == 1)).sum())}/'
                                f'{int(((m.recommendation_ref == 1) & (m.recommendation == 0)).sum())}')

row["causalpfn_seconds"] = cfg.get("causalpfn_seconds")

new = pd.DataFrame([row])
if SUMMARY.exists():
    old = pd.read_csv(SUMMARY)
    old = old[old.run_label != label]
    out = pd.concat([old, new], ignore_index=True)
else:
    out = new
out.to_csv(SUMMARY, index=False)
print(new.T.to_string(header=False))
