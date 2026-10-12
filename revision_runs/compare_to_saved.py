"""Compare a revision run's est_effect with the saved paper export.

Usage: python compare_to_saved.py RUN_CSV [SAVED_CSV]
Default SAVED_CSV is BRACE/Datasets/metabric_evaluation_data_20260617.csv
(594 METABRIC test patients, est_effect from the paper's seed 123 PAM50 run).
Writes RUN_CSV with suffix _vs_saved.json next to it.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SAVED_DEFAULT = Path("/Users/sindypinero/Downloads/4_Deanna/BRACE/Datasets/metabric_evaluation_data_20260617.csv")


def main():
    run_csv = Path(sys.argv[1])
    saved_csv = Path(sys.argv[2]) if len(sys.argv) > 2 else SAVED_DEFAULT

    run = pd.read_csv(run_csv)
    saved = pd.read_csv(saved_csv, usecols=["patientID", "time.x", "event.x", "treatment.x", "est_effect", "recommendation"])
    saved = saved.rename(columns={"time.x": "time", "event.x": "event", "treatment.x": "treatment"})

    ids_run, ids_saved = set(run.patientID), set(saved.patientID)
    m = saved.merge(run, on="patientID", suffixes=("_saved", "_new"))

    out = {
        "run_csv": str(run_csv),
        "saved_csv": str(saved_csv),
        "n_run": len(run),
        "n_saved": len(saved),
        "n_common_patients": len(m),
        "n_only_in_run": len(ids_run - ids_saved),
        "n_only_in_saved": len(ids_saved - ids_run),
        "same_row_order": bool(len(run) == len(saved) and (run.patientID.values == saved.patientID.values).all()),
    }
    if len(m) >= 3:
        for col in ["time", "event", "treatment"]:
            out[f"{col}_mismatches"] = int((~np.isclose(m[f"{col}_saved"], m[f"{col}_new"])).sum())
        a, b = m.est_effect_saved.to_numpy(), m.est_effect_new.to_numpy()
        out.update({
            "pearson": float(stats.pearsonr(a, b)[0]),
            "spearman": float(stats.spearmanr(a, b)[0]),
            "max_abs_diff": float(np.max(np.abs(a - b))),
            "mean_abs_diff": float(np.mean(np.abs(a - b))),
            "saved_est_effect_mean_sd": [float(a.mean()), float(a.std(ddof=1))],
            "new_est_effect_mean_sd": [float(b.mean()), float(b.std(ddof=1))],
            "sign_agreement": float(np.mean(np.sign(a) == np.sign(b))),
            "split_saved_no_chemo_chemo": [int((m.recommendation_saved == 0).sum()), int((m.recommendation_saved == 1).sum())],
            "split_new_no_chemo_chemo": [int((m.recommendation_new == 0).sum()), int((m.recommendation_new == 1).sum())],
            "recommendation_crosstab_saved_rows_new_cols": pd.crosstab(m.recommendation_saved, m.recommendation_new).to_dict(),
        })
        out["reproduces"] = bool(
            out["n_common_patients"] == 594 and out["pearson"] >= 0.99 and out["spearman"] >= 0.99
            and out["split_new_no_chemo_chemo"] == [524, 70]
        )
    else:
        out["reproduces"] = False

    out_path = run_csv.with_name(run_csv.stem + "_vs_saved.json")
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
