"""Run E (22 inputs x 7 seeds) and G (7 training subsamples of PAM50 clean + clinical, seed 123 split).

Usage: python run_EG.py [N_WORKERS]
Resumable: a run is skipped when LABEL.json, LABEL_eval.json and LABEL_policy_value.json all exist.
Each run: master.R -> evaluate_run.R -> policy_value_for_run.py (B = 500). Failures are logged and the
remaining runs continue. Afterwards the summary rows are appended in a fixed order (append_summary.py).
Logs: revision_runs/logs/E/LABEL.log, revision_runs/logs/G/LABEL.log, revision_runs/logs/run_EG.log
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from brace_paths import PYTHON as PY, REV

os.environ["BRACE_PYTHON"] = PY  # R child processes call CausalPFN with the same interpreter
CLEAN = REV / "7_PAM50_clean.csv"
SEEDS = [123, 26, 314, 11, 24, 150, 111]  # Supplementary Table S2
B_POLICY = 500

# 21 biological inputs of Table 3 (PAM50 with the clean 50 gene list) plus PAM50 clean + clinical (C3 spec)
INPUTS = [
    ("pearson100", "Pearson corr. coeff. (100)", "FS_CHOICE <- 'correlation'"),
    ("gs1_driver_core", "Driver Core", "GENE_SET <- '1'"),
    ("gs2_driver_ext_wgs", "Driver Extended WGS", "GENE_SET <- '2'"),
    ("gs3_driver_found", "Driver Foundation", "GENE_SET <- '3'"),
    ("gs4A_tcga_cohorts", "TCGA BRCA Cohorts", "GENE_SET <- '4A'"),
    ("gs4B_tcga_subgroups", "TCGA BRCA Subgroups", "GENE_SET <- '4B'"),
    ("gs4C_tcga_all", "TCGA BRCA All", "GENE_SET <- '4C'"),
    ("gs5_interval_screen", "Interval Screen", "GENE_SET <- '5'"),
    ("gs6_prolif_core_cc", "Prolif Core CC", "GENE_SET <- '6'"),
    ("gs7_pam50_clean", "PAM50 (clean, 50 genes)", f"GENE_FILE <- '{CLEAN}'"),
    ("gs8_mammaprint", "MammaPrint", "GENE_SET <- '8'"),
    ("gs9_blueprint", "Blueprint", "GENE_SET <- '9'"),
    ("gs10_oncotype", "Oncotype", "GENE_SET <- '10'"),
    ("gs11_endopredict", "EndoPredict", "GENE_SET <- '11'"),
    ("gs12_index_bci", "Index BCI", "GENE_SET <- '12'"),
    ("ct_I_k", "Control theory Type I, K", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'I'; CT_DIRECTION <- ''"),
    ("ct_I_kin", "Control theory Type I, K_in", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'I'; CT_DIRECTION <- 'in'"),
    ("ct_I_kout", "Control theory Type I, K_out", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'I'; CT_DIRECTION <- 'out'"),
    ("ct_II_k", "Control theory Type II, K", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'II'; CT_DIRECTION <- ''"),
    ("ct_II_kin", "Control theory Type II, K_in", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'II'; CT_DIRECTION <- 'in'"),
    ("ct_II_kout", "Control theory Type II, K_out", "FS_CHOICE <- 'control_theory'; CT_TYPE <- 'II'; CT_DIRECTION <- 'out'"),
    ("pam50clean_clinical", "PAM50 clean + clinical (C3 spec)", f"GENE_FILE <- '{CLEAN}'; USE_CLINICAL <- TRUE"),
]


def jobs():
    out = []
    for key, desc, settings in INPUTS:
        for s in SEEDS:
            out.append(dict(label=f"E_{key}_s{s}", group="E", desc=desc, ref="",
                            settings=f"SEED <- {s}; {settings}"))
    for s in SEEDS:
        out.append(dict(label=f"G_pam50clean_clinical_sub80_s{s}", group="G",
                        desc="PAM50 clean + clinical, 80% training subsample", ref="C3_pam50_clinical_seed123",
                        settings=f"SEED <- 123; TRAIN_SUBSAMPLE <- 0.8; SUBSAMPLE_SEED <- {s}; "
                                 f"GENE_FILE <- '{CLEAN}'; USE_CLINICAL <- TRUE"))
    return out


def done(label):
    return all((REV / f"{label}{suf}").exists() for suf in (".json", "_eval.json", "_policy_value.json"))


def run_one(job):
    label = job["label"]
    log = REV / "logs" / job["group"] / f"{label}.log"
    t0 = time.time()
    with open(log, "w") as fh:
        steps = [
            ["Rscript", "-e", f"RUN_LABEL <- '{label}'; {job['settings']}; source('master.R')"],
            ["Rscript", "evaluate_run.R", label],
            [PY, str(REV / "policy_value_for_run.py"), label, str(B_POLICY)],
        ]
        for cmd in steps:
            fh.write(f"\n$ {' '.join(cmd)}\n"); fh.flush()
            r = subprocess.run(cmd, cwd=REV / "R", stdout=fh, stderr=subprocess.STDOUT)
            if r.returncode != 0:
                return label, False, time.time() - t0
    return label, True, time.time() - t0


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    for g in ("E", "G"):
        (REV / "logs" / g).mkdir(parents=True, exist_ok=True)
    all_jobs = jobs()
    todo = [j for j in all_jobs if not done(j["label"])]
    mlog = open(REV / "logs" / "run_EG.log", "a")
    say = lambda m: (mlog.write(time.strftime("%H:%M:%S ") + m + "\n"), mlog.flush(), print(m, flush=True))
    say(f"{len(all_jobs)} runs, {len(all_jobs) - len(todo)} already done, {len(todo)} to run, {workers} workers")
    failed = []
    with ThreadPoolExecutor(workers) as ex:
        futs = [ex.submit(run_one, j) for j in todo]
        for i, f in enumerate(as_completed(futs), 1):
            label, ok, secs = f.result()
            say(f"[{i}/{len(todo)}] {label} {'ok' if ok else 'FAILED'} ({secs:.0f} s)")
            if not ok:
                failed.append(label)
    for j in all_jobs:  # summary rows in fixed order
        if done(j["label"]):
            args = [PY, str(REV / "append_summary.py"), j["label"], j["desc"]] + ([j["ref"]] if j["ref"] else [])
            subprocess.run(args, stdout=subprocess.DEVNULL, check=True)
    say(f"finished; failed: {failed if failed else 'none'}")


if __name__ == "__main__":
    main()
