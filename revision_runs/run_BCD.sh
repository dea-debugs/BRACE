#!/bin/zsh
# Runs B, C (3 variants) and D (2 outcomes) in sequence; stops at the first failure.
# Each run: master.R (CausalPFN, writes LABEL.csv/.json) -> evaluate_run.R (stratified RMST at 150)
#           -> policy_value_for_run.py (overlap policy value, B = 500) -> append_summary.py
set -e
REV=${0:A:h}                       # this folder (revision_runs); data paths come from BRACE_ROOT, see R/paths.R
PY=${BRACE_PYTHON:-$(command -v python)}
export BRACE_PYTHON=$PY            # R calls CausalPFN with the same interpreter
CLEAN=$REV/7_PAM50_clean.csv
cd $REV/R

run() {  # label, R settings, input description, reference label
  local label=$1 settings=$2 desc=$3 ref=$4
  echo "=== $label ($(date +%H:%M:%S))"
  Rscript -e "RUN_LABEL <- '$label'; SEED <- 123; $settings; source('master.R')" > $REV/$label.log 2>&1 \
    || { echo "master.R failed for $label, see $label.log"; tail -20 $REV/$label.log; exit 1; }
  grep "CausalPFN call" $REV/$label.log
  Rscript evaluate_run.R $label 2>/dev/null
  $PY $REV/policy_value_for_run.py $label 500 > $REV/${label}_policy_value.log 2>&1 \
    || { echo "policy value failed for $label"; tail -20 $REV/${label}_policy_value.log; exit 1; }
  $PY $REV/append_summary.py $label "$desc" $ref > /dev/null
}

run B_pam50clean_seed123          "GENE_FILE <- '$CLEAN'"                                  "PAM50 clean list (50 genes)"           A_baseline_pam50_seed123
run C1_clinical_only_seed123      "GENE_SET <- 'none'; USE_CLINICAL <- TRUE"               "clinical only (10 + missing ind.)"     A_baseline_pam50_seed123
run C2_pam50_only_seed123         "GENE_FILE <- '$CLEAN'"                                  "PAM50 only (50 genes)"                 A_baseline_pam50_seed123
run C3_pam50_clinical_seed123     "GENE_FILE <- '$CLEAN'; USE_CLINICAL <- TRUE"            "PAM50 (50) + clinical"                 A_baseline_pam50_seed123
run D1_pam50_clinical_ipcw150_seed123   "GENE_FILE <- '$CLEAN'; USE_CLINICAL <- TRUE; OUTCOME <- 'ipcw150'"        "PAM50 (50) + clinical" C3_pam50_clinical_seed123
run D2_pam50_clinical_pseudo150_seed123 "GENE_FILE <- '$CLEAN'; USE_CLINICAL <- TRUE; OUTCOME <- 'pseudo_rmst150'" "PAM50 (50) + clinical" C3_pam50_clinical_seed123
echo "=== all done ($(date +%H:%M:%S))"
