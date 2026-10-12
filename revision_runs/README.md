# BRACE revision runs (JBI)

This folder holds the code, the settings and the summary outputs behind the revised BRACE paper. `REPORT.md` describes every run and its results in detail. This README explains how to rebuild the software environment, how the pieces connect, and how to reproduce each run, table and figure.

## What is here and what is not

Included: all Python and shell scripts, the edited R pipeline in `R/`, the clean PAM50 list, the clinical covariate table, the summary tables of runs A to G, the run F summary, the figures and the report. The policy value analysis that the revision reuses sits in `../policy_value_analysis/`.

Not included: the METABRIC expression data, the cBioPortal clinical files, the R data cache, the per run exchange folders, logs, and the per run CSV, JSON and PDF outputs. All of these can be rebuilt with the steps below.

## 1. Recreate the conda environment

The runs used a native Apple silicon (osx arm64) conda environment called `brace` with Python 3.12.15. `environment.txt` lists `pip show causalpfn` followed by the full `pip freeze`.

```bash
CONDA_SUBDIR=osx-arm64 conda create -n brace python=3.12 -y
conda activate brace
conda config --env --set subdir osx-arm64   # keeps later installs native
# install the exact versions recorded in environment.txt
sed -n '/=== pip freeze ===/,$p' environment.txt | grep '==' > requirements_brace.txt
pip install -r requirements_brace.txt
```

On other platforms, drop the `CONDA_SUBDIR` and `subdir` lines. The key versions are causalpfn 0.1.4, torch 2.14.1, numpy 2.5.4, pandas, scikit-learn and lifelines.

### CausalPFN checkpoint

On first use CausalPFN downloads its weights from Hugging Face (`vdblm/causalpfn`). The runs used:

| File | Size | SHA256 |
|---|---|---|
| `causalpfn_v0.pt` | 75,419,490 bytes | `4f0f53717f1c62426fe98faa0459d33bc8f9862a87e8484ca3869ca899e7841d` |

Check your copy with `shasum -a 256 causalpfn_v0.pt` in the Hugging Face cache folder.

### R

The R side used R 4.5.3 with survival, survRM2, grf 2.6.1, caret, zeallot, FSelectorRcpp, randomForestSRC, jsonlite and causalTree (installed from GitHub). With R 4.5 and recent compilers causalTree needs extra flags. Use `R/Makevars.causalTree` through `R_MAKEVARS_USER` when installing it.

### Paths

The scripts use absolute paths from the machine where the runs were made (`/Users/sindypinero/Downloads/4_Deanna/...`). Before running anything elsewhere, edit these constants:

| File | Constants |
|---|---|
| `R/hte.R` | `CAUSALPFN_PYTHON`, `CAUSALPFN_SCRIPT`, `CAUSALPFN_IO_DIR` |
| `R/master.R` | `CLINICAL_FILE`, `REVISION_DIR`, `DATA_DIR`, `GENE_SET_DIR` |
| `R/feature_selection.R` | default `path` arguments |
| `R/evaluate_run.R` | `REVISION_DIR` |
| `R/revision_helpers.R` | path to `censoring_outcomes.py` |
| Python scripts and `run_BCD.sh` | `REV`, `BRACE`, `BASE` and `PY` near the top |

The expression data are read from `DATA_DIR/dataset_genes.csv`, which is built by `R/1_create_nice_dataset.R`. Gene sets are read from `GENE_SET_DIR`.

## 2. How R calls CausalPFN

`R/hte.R` writes two files into a per run folder `causalpfn_io/RUN_LABEL/`:

* `causalpfn_train.csv`: columns `Y` (outcome), `W` (treatment, 0 or 1), then the features
* `causalpfn_test.csv`: the same feature columns

It then runs

```
cd causalpfn_io/RUN_LABEL && <CAUSALPFN_PYTHON> <CAUSALPFN_SCRIPT>
```

`run_causalpfn.py` checks that the feature columns match, fits `CATEEstimator` with default settings (seeds 123 for random, numpy and torch), and writes `causalpfn_predictions.csv` with one column, `risk_score`, the estimated effect for each test patient. R reads that file back as `est_effect` and stops if the call fails. The script sets `KMP_DUPLICATE_LIB_OK=TRUE` and `OMP_NUM_THREADS=1` before importing anything, and it ends with `os._exit(0)`. Without these, torch, faiss and scikit learn load clashing OpenMP libraries, and the process can hang or crash on exit. CausalPFN is deterministic for fewer than 4,096 rows, so the seed of a run only changes the train and test split.

## 3. Reproducing the runs

Every run follows the same chain:

1. `R/master.R` builds the split, selects features, calls CausalPFN, and writes `LABEL.csv` and `LABEL.json`.
2. `Rscript R/evaluate_run.R LABEL` computes the RMST difference (Followed minus Not Followed) at 150 months, pooled and within each observed arm, and writes `LABEL_eval.json`.
3. `python policy_value_for_run.py LABEL 500` gives the doubly robust overlap policy values with 500 bootstrap replicates and writes `LABEL_policy_value.json`.
4. `python append_summary.py LABEL "description" REFERENCE_LABEL` adds the run as a row of `summary.csv`.

Settings are R variables set before `source("master.R")`: `RUN_LABEL`, `SEED`, `GENE_SET`, `GENE_FILE`, `OUTCOME`, `USE_CLINICAL`, `FS_CHOICE`, `CT_TYPE`, `CT_DIRECTION`, `TRAIN_SUBSAMPLE` and `SUBSAMPLE_SEED`. The defaults reproduce run A. All commands below run from `R/` with the `brace` Python as `$PY`.

### Inputs prepared once

* Clinical files: `data_clinical_sample.txt` and `data_clinical_patient.txt` from the cBioPortal datahub (`public/brca_metabric`), saved in `cbioportal/`. The SHA256 sums used were `9a2291e065e2a0a5337b4eaeb09836e70b1bc3e1ed17714a013b53515697e6dc` (sample) and `2a6bbadc6e4b2b07c0505559d438da959fc5c9aba55c7194abc394e2d6454995` (patient).
* `python build_clinical_covariates.py` writes `clinical_covariates.csv` with 10 raw clinical covariates. Median imputation with missing indicators is fitted on the training partition inside `master.R`.
* `policy_value_for_run.py` executes the data section of `../policy_value_analysis/policy_value.py` unchanged, apart from its two input paths. It also needs the original analysis extracts in `BRACE/Datasets/analysis_extracts/`.

### Run A: baseline reproduction

```bash
Rscript master.R                                   # defaults: original PAM50 list (GENE_SET 7), seed 123
Rscript evaluate_run.R A_baseline_pam50_seed123
$PY ../policy_value_for_run.py A_baseline_pam50_seed123 500
$PY ../compare_to_saved.py ../A_baseline_pam50_seed123.csv   # compares with the paper export
```

### Runs B, C and D

`zsh ../run_BCD.sh` runs them in sequence and stops at the first failure:

* B: clean 50 gene PAM50 list
* C1: clinical only
* C2: PAM50 only
* C3: PAM50 plus clinical, the primary specification of the revision
* D1: C3 with the IPCW outcome at 150 months
* D2: C3 with the jackknife pseudo observation of RMST at 150 months

The D outcomes are computed on the training partition by `censoring_outcomes.py`.

### Runs E and G

`$PY ../run_EG.py 5` runs everything with 5 parallel workers and can be resumed after an interruption:

* E: 22 inputs × 7 seeds (123, 26, 314, 11, 24, 150, 111), 154 runs
* G: seven refits of the C3 specification, each on an 80% subsample of the seed 123 training partition, with the test set fixed

Then:

```bash
$PY ../summarise_E.py        # E_seed_summary.csv, E_seed_long.csv, E_pam50_paired.csv
$PY ../rank_E_policy.py      # E_policy_ranking.csv
$PY ../analyse_G.py          # G_patient_spread.csv, G_summary.json, G_confident_consensus.csv
for l in G_confident_consensus C3_restricted_to_G_confident; do
  Rscript evaluate_run.R $l
  $PY ../policy_value_for_run.py $l 500 ../C3_pam50_clinical_seed123.csv
done
```

### Run F: semi synthetic benchmark

```bash
$PY ../F_semisynthetic.py --smoke          # one setting, one replicate, to time it
nohup $PY ../F_semisynthetic.py --threads 4 > /dev/null 2>&1 &
```

The full grid has 36 settings with 10 replicates each and took a little over one hour on an Apple silicon Mac. It can be resumed: replicates already on disk are skipped. Progress goes to `logs/F/F_run.log`, and the result is `F_semisynthetic/F_summary.csv`. The causal survival forest comparator is `R/csf_fit.R`.

### Figures

`$PY make_figures.py` writes both PNG files to `figures/` at 300 dpi.

## 4. Where each table and figure comes from

| Revised paper item | Content | Script | Output file |
|---|---|---|---|
| Methods, software | Environment and checkpoint | | `environment.txt` |
| Baseline reproduction (run A) | Paper results reproduced from code | `R/master.R`, `R/evaluate_run.R`, `compare_to_saved.py` | `summary.csv` row `A_baseline_pam50_seed123` |
| UBE2T correction (run B) | Effect of the clean PAM50 list | `run_BCD.sh` | `summary.csv` row `B_pam50clean_seed123` |
| Clinical covariates table (run C) | Clinical only, PAM50 only, PAM50 plus clinical | `build_clinical_covariates.py`, `run_BCD.sh` | `summary.csv` rows `C1`, `C2`, `C3` |
| Figure: Kaplan Meier curves at 150 months | Primary specification C3: pooled, no chemotherapy arm, chemotherapy arm | `make_figures.py` | `figures/C3_km_150.png` |
| Censoring aware outcomes table (run D) | IPCW and pseudo observation outcomes | `censoring_outcomes.py`, `run_BCD.sh` | `summary.csv` rows `D1`, `D2` |
| Policy value results | Overlap policy value versus observed care and treat none | `policy_value_for_run.py`, `../policy_value_analysis/policy_value.py` | `LABEL_policy_value.json` per run, collected in `summary.csv`; original analysis in `../policy_value_analysis/policy_value_results.json` |
| Table 3 revised with seed ranges (run E) | Pooled RMST difference across seven splits for every input | `run_EG.py`, `summarise_E.py` | `E_seed_summary.csv`, `E_seed_long.csv`, `E_pam50_paired.csv` |
| Policy value ranking table (run E) | Inputs ranked on V(BRACE) minus V(observed care), with reference policies | `rank_E_policy.py` | `E_policy_ranking.csv` |
| Figure: policy value ranking | Dot and range plot of the ranking | `make_figures.py` | `figures/E_policy_ranking.png` |
| Supplementary Table S2 | Seeds and arm counts for each split | `run_EG.py`, `summarise_E.py` | `E_seed_long.csv` (columns `n_test_chemo`, `n_test_no_chemo`) |
| Recommendation intervals and abstention (run G) | Per patient spread across seven refits, confident subset | `run_EG.py`, `analyse_G.py` | `G_patient_spread.csv`, `G_summary.json`, `summary.csv` row `G_confident_consensus` |
| Semi synthetic benchmark table (run F) | PEHE, sign agreement and regret for CausalPFN and the causal survival forest | `F_semisynthetic.py`, `R/csf_fit.R` | `F_semisynthetic/F_summary.csv` |

Table numbers in the revised manuscript may change. The content column says which result each row covers, and `REPORT.md` has the matching numbers section by section.
