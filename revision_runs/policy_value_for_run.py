"""Overlap population policy value for one revision run, reusing policy_value_analysis/policy_value.py.

Usage: python policy_value_for_run.py RUN_LABEL [B] [UNIVERSE_CSV]
UNIVERSE_CSV (optional) is a run CSV whose patients form the full test partition; nuisance models are fitted
on all other patients. Default: the run's own patients. For the seed 123 partition the original test/train
extracts are used unchanged; for other partitions (other seeds) test and train rows come from the
all1980 extract of the same script. Restricting RUN_LABEL.csv to a subset (run G) restricts evaluation only.
The data/helper section of policy_value.py (everything before "# ---------------- 0. validation") is executed
unchanged, except that its two input paths (U, S) are pointed at this Mac. The run's recommendations replace
te['R']; then run('clinical+pam50', 'train', overlap=(0.1, 0.9)) gives the overlap population policy values,
exactly as in the original analysis, and the original bootstrap scheme (test and training resampled) gives
percentile CIs. Writes revision_runs/RUN_LABEL_policy_value.json.
"""
import json
import sys

import numpy as np
import pandas as pd

from brace_paths import DATASETS_DIR, PV_SCRIPT, REV

label = sys.argv[1]
B = int(sys.argv[2]) if len(sys.argv) > 2 else 200
universe_csv = sys.argv[3] if len(sys.argv) > 3 else None

# policy_value.py reads S + 'patient.txt' and S + 'sample.txt'; provide them under those names
pv_inputs = REV / "policy_value_inputs"
pv_inputs.mkdir(exist_ok=True)
for name, src in (("patient.txt", "data_clinical_patient.txt"), ("sample.txt", "data_clinical_sample.txt")):
    dst = pv_inputs / name
    if not dst.exists():
        dst.symlink_to(REV / "cbioportal" / src)

src = PV_SCRIPT.read_text()
head = src.split("# ---------------- 0. validation")[0]
old_U = "U = '/mnt/user-data/uploads/4_Deanna/BRACE/Datasets/analysis_extracts/'"
old_S = "S = '/tmp/claude-0/-home-claude/b1f82c5d-573c-5bd8-b091-63beff174c30/scratchpad/pv/'"
assert old_U in head and old_S in head
head = head.replace(old_U, f"U = '{DATASETS_DIR}/analysis_extracts/'").replace(old_S, f"S = '{pv_inputs}/'")
head = head.replace("OUT = sys.argv[1] if len(sys.argv) > 1 else S", "OUT = S")
ns = {"__name__": "policy_value"}
exec(compile(head, str(PV_SCRIPT), "exec"), ns)

te, tr, al, run = ns["te"], ns["tr"], ns["al"], ns["run"]
rec = pd.read_csv(REV / f"{label}.csv")[["patientID", "recommendation"]]
universe = set(pd.read_csv(universe_csv).patientID) if universe_csv else set(rec.patientID)
assert set(rec.patientID) <= universe
if universe == set(te.patientID):  # seed 123 partition: original extracts
    te_src, tr = te.drop(columns=["R", "recommendation", "followed"]), tr
    partition = "seed 123 (original test594 / train1386 extracts)"
else:
    te_src = al[al.patientID.isin(universe)].copy()
    tr = al[~al.patientID.isin(universe)].reset_index(drop=True)
    partition = "from all1980 extract (run's own test partition)"
assert len(tr) == 1386, len(tr)
te_new = te_src.merge(rec, on="patientID", how="inner", validate="one_to_one").reset_index(drop=True)
assert len(te_new) == len(rec)
te_new["R"] = te_new.recommendation.astype(int)
ns["te"], ns["tr"] = te_new, tr  # run() defaults to the module level te/tr

OV = (0.1, 0.9)
pv_ov, _ = run("clinical+pam50", "train", te_=te_new, tr_=tr, overlap=OV)
pv_full, _ = run("clinical+pam50", "train", te_=te_new, tr_=tr)

boot = {k: [] for k in ["ov_BRACE_DR", "ov_BRACE_IPW", "ov_obs", "ov_none_DR", "ov_none_IPW", "ov_all_DR"]}
for b in range(B):
    teb = te_new.sample(len(te_new), replace=True, random_state=b).reset_index(drop=True)
    trb = tr.sample(len(tr), replace=True, random_state=10_000 + b).reset_index(drop=True)
    try:
        pvo, _ = run("clinical+pam50", "train", te_=teb, tr_=trb, overlap=OV)
    except Exception:
        continue
    boot["ov_BRACE_DR"].append(pvo["BRACE"]["DR"]); boot["ov_BRACE_IPW"].append(pvo["BRACE"]["IPW"])
    boot["ov_obs"].append(pvo["observed"]["IPW"]); boot["ov_none_DR"].append(pvo["treat_none"]["DR"])
    boot["ov_none_IPW"].append(pvo["treat_none"]["IPW"]); boot["ov_all_DR"].append(pvo["treat_all"]["DR"])
bt = {k: np.array(v) for k, v in boot.items()}


def ci(a):
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


out = {
    "run_label": label,
    "adjustment_set": "clinical+pam50 (policy_value.py SETS, 49 PAM50 genes as in the original analysis)",
    "nuisance": "fitted on training partition",
    "partition": partition,
    "n_evaluated": len(te_new),
    "horizon": float(ns["H"]),
    "overlap": list(OV),
    "n_overlap": pv_ov["n_overlap"],
    "reclassified_in_overlap": pv_ov["reclassified_in_overlap"],
    "overlap_policy_values": {k: v for k, v in pv_ov.items() if isinstance(v, dict)},
    "overlap_contrasts": {
        "V(BRACE) - V(observed) DR": pv_ov["BRACE"]["DR"] - pv_ov["observed"]["IPW"],
        "V(BRACE) - V(treat none) DR": pv_ov["BRACE"]["DR"] - pv_ov["treat_none"]["DR"],
        "V(BRACE) - V(treat none) IPW": pv_ov["BRACE"]["IPW"] - pv_ov["treat_none"]["IPW"],
    },
    "full_population_policy_values": pv_full,
    "bootstrap_B": len(bt["ov_obs"]),
    "bootstrap_CI": {
        "overlap: V(BRACE) DR": ci(bt["ov_BRACE_DR"]),
        "overlap: V(BRACE) - V(observed) DR": ci(bt["ov_BRACE_DR"] - bt["ov_obs"]),
        "overlap: V(BRACE) - V(treat none) DR": ci(bt["ov_BRACE_DR"] - bt["ov_none_DR"]),
        "overlap: V(BRACE) - V(treat none) IPW": ci(bt["ov_BRACE_IPW"] - bt["ov_none_IPW"]),
    },
}
(REV / f"{label}_policy_value.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps({k: out[k] for k in ["n_overlap", "overlap_contrasts", "bootstrap_CI"]}, indent=1, default=float))
print("V(BRACE) DR overlap:", out["overlap_policy_values"]["BRACE"]["DR"])
