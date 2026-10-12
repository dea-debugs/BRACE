# BRACE revision runs: report

Prepared 11 October 2026 on the Mac (Apple M4 Pro). All files are in `revision_runs/`. The original exports, R code and gene sets were not modified; all R edits are in copies under `revision_runs/R/`.

Conventions used throughout:

* RMST difference means Followed minus Not Followed at 150 months, in months, computed with survRM2 (the estimator of Uno et al. used in the paper). "No chemo arm" and "chemo arm" are the treatment stratified comparisons within each observed arm. "Pooled" is the headline comparison over all test patients.
* Policy value means the restricted mean survival time at 150 months that a recommendation policy would achieve, estimated with the doubly robust (DR) estimator of `policy_value_analysis/policy_value.py` on the overlap population (patients with estimated propensity between 0.1 and 0.9), with nuisance models fitted on the training partition. Intervals are percentile intervals from 500 bootstrap draws, using the same scheme as the original analysis.
* Every run writes `LABEL.csv` (patientID, treatment, time, event, est_effect, recommendation, run_label), `LABEL.json` (settings and seed), `LABEL_eval.json` and `LABEL_policy_value.json`. One row per run is in `summary.csv` (169 rows).

---

## 1. Setup

| Item | Result |
|---|---|
| Environment | conda env `brace`, native arm64, Python 3.12.15, `subdir: osx-arm64` pinned in the env |
| CausalPFN | causalpfn 0.1.4, torch 2.14.1, numpy 2.5.4 (`environment.txt` has pip show and pip freeze) |
| Checkpoint | `causalpfn_v0.pt` from Hugging Face `vdblm/causalpfn`, 75,419,490 bytes, SHA256 `4f0f53717f1c62426fe98faa0459d33bc8f9862a87e8484ca3869ca899e7841d` |
| Bridge script | `run_causalpfn.py` recreated from the interface in `R Code/hte.R` (reads `causalpfn_train.csv` with Y, W, features and `causalpfn_test.csv`; writes `causalpfn_predictions.csv` with `risk_score`); CATEEstimator defaults, seed 123 |
| PAM50 list | `7_PAM50_clean.csv`: LF line endings, 50 genes, stray byte 0xCA removed from UBE2T; original untouched |
| R | R 4.5.3 with grf 2.6.1, survRM2, caret, zeallot, FSelectorRcpp, randomForestSRC and causalTree (GitHub; build recipe in `R/Makevars.causalTree`) |
| Speed | expression data cached as `cache/dataset_genes.rds`; one full run takes about 8 seconds, of which CausalPFN takes 3 to 5 seconds |

The first environment was an Intel (x86_64) Python under Rosetta, which could only install an old torch build that cannot work with NumPy 2, so the environment was rebuilt as native arm64. On arm64 the three OpenMP runtimes bundled with torch, faiss and scikit learn crash when they run threads together, so `run_causalpfn.py` limits OpenMP to one thread and ends with `os._exit(0)`; it runs in seconds and returns control to R reliably. The bridge script is not yet committed to the BRACE GitHub repository, which is still needed to fully close item 24.

Checklist items closed: 24 (missing `run_causalpfn.py` and environment record, pending the GitHub commit) and the checkpoint part of 34.

---

## 2. Run A: baseline reproduction

| Check | Paper export | This run |
|---|---|---|
| Test patients | 594 | same 594 patients, same order |
| Pearson / Spearman correlation of est_effect | | 0.99999999999790 / 0.9999999999999999 |
| Largest absolute difference in est_effect | | 0.00014 |
| Recommendations (no chemo / chemo) | 524 / 70 | 524 / 70, no patient changes |
| Pooled RMST difference | 18.94 (9.31 to 28.57) | 18.94 (9.31 to 28.57), p = 0.0001 |
| No chemo arm | 19.19, p = 0.018 | 19.19 (3.31 to 35.08), p = 0.018 |
| Chemo arm | 6.24, p = 0.58 | 6.24 (−16.11 to 28.59), p = 0.58 |
| Overlap policy value, V(BRACE) minus observed | 0.50 | 0.50 (−17.03 to 10.64), n = 173 |

The published pipeline reproduces exactly: same split, same estimates up to float rounding, same recommendations and the same survival and policy value results. This confirms that the recreated `run_causalpfn.py`, the checkpoint and the R pipeline are the ones used for the paper. Every later run differs from A only in the setting named in its section.

Checklist items closed: none directly; this is the prerequisite for all others.

---

## 3. Run B: UBE2T corrected

| | A (49 genes) | B (50 genes, clean list) |
|---|---|---|
| Recommended chemo | 70 | 69 |
| Recommendations changed versus A | | 7 of 594 (3 to chemo, 4 to no chemo) |
| No chemo arm | 19.19 (3.31 to 35.08), p = 0.018 | 13.73 (−1.84 to 29.31), p = 0.084 |
| Chemo arm | 6.24 (−16.11 to 28.59), p = 0.58 | 6.76 (−15.56 to 29.08), p = 0.55 |
| Pooled | 18.94 (9.31 to 28.57), p = 0.0001 | 17.27 (7.64 to 26.89), p = 0.0004 |
| V(BRACE) minus observed | 0.50 (−17.03 to 10.64) | 0.94 (−18.03 to 13.62) |

Adding the missing gene changes only 7 recommendations, but the stratified result in the no chemo arm, which the paper reports as significant, is no longer significant. The pooled difference stays clearly positive. The headline numbers in the paper therefore depend on the file error, and the corrected list should be used for any reported PAM50 result.

Checklist items closed: 34 (UBE2T).

---

## 4. Run C: clinical covariates

| Input (seed 123) | Features | Chemo recs | Changed vs A | No chemo arm | Chemo arm | Pooled | V(BRACE) minus observed | V(BRACE) minus treat none |
|---|---|---|---|---|---|---|---|---|
| Clinical only | 15 | 39 | 93 | 46.86 (19.64 to 74.09), p = 0.001 | −10.12 (−34.30 to 14.06) | 17.36 (6.82 to 27.90) | −2.96 (−16.79 to 11.91) | 3.22 (−5.76 to 13.44) |
| PAM50 only (same as B) | 50 | 69 | 7 | 13.73 (−1.84 to 29.31), p = 0.084 | 6.76 (−15.56 to 29.08) | 17.27 (7.64 to 26.89) | 0.94 (−18.03 to 13.62) | 7.12 (−11.79 to 17.21) |
| PAM50 plus clinical | 65 | 105 | 91 | 17.89 (4.64 to 31.14), p = 0.008 | −12.34 (−31.89 to 7.21) | 13.58 (4.12 to 23.03) | −0.40 (−18.62 to 14.61) | 5.78 (−11.99 to 17.92) |

Clinical covariates are age, positive lymph nodes, NPI, ER (IHC), menopausal state, grade, tumour size, tumour stage, HER2 and PR (the last five from cBioPortal `data_clinical_sample.txt`, joined on PATIENT_ID), with median imputation and missing indicators fitted on the training partition only (15 columns). Adding them to PAM50 changes 91 recommendations, keeps a significant separation in the no chemo arm, and lowers the pooled difference. The clinical only result in the no chemo arm looks large, but only 11 untreated patients are in its Not Followed group, and no variant has a policy value interval that excludes zero.

Checklist items closed: 1 (clinical covariates in X).

---

## 5. Run D: censoring aware outcomes

| Outcome supplied to CausalPFN (PAM50 plus clinical) | Chemo recs | Changed vs C | No chemo arm | Chemo arm | Pooled | V(BRACE) minus observed | V(BRACE) minus treat none |
|---|---|---|---|---|---|---|---|
| min(U, 360), as in C | 105 | | 17.89 (4.64 to 31.14) | −12.34 (−31.89 to 7.21) | 13.58 (4.12 to 23.03) | −0.40 (−18.62 to 14.61) | 5.78 (−11.99 to 17.92) |
| IPCW pseudo outcome at 150 months | 255 | 178 | 2.29 (−6.34 to 10.92) | −13.13 (−31.75 to 5.49) | −1.29 (−9.27 to 6.69) | −12.78 (−38.26 to 9.08) | −6.60 (−33.25 to 16.67) |
| Jackknife pseudo observation of RMST at 150 months | 279 | 190 | 2.16 (−5.99 to 10.31) | 7.36 (−11.20 to 25.91) | 4.18 (−3.49 to 11.85) | 3.48 (−17.09 to 22.78) | 9.66 (−13.47 to 26.87) |

The IPCW outcome uses a Kaplan Meier censoring curve fitted by treatment arm on the training partition (273 of 1,386 training patients get zero weight; the smallest censoring survival was 0.48, so no weights were clipped), and the pseudo observations are computed in Python on the whole training partition. With either censoring aware outcome the model recommends chemotherapy to far more patients and the survival separation disappears in both arms and pooled. This suggests that the separation reported with min(U, 360) is tied to that outcome encoding rather than being a robust property of the recommendations.

Checklist items closed: 2 (censoring aware outcomes).

---

## 6. Run E: seven seeds for every input

Seeds 123, 26, 314, 11, 24, 150 and 111 from Supplementary Table S2; each seed gives a different 70/30 split. PAM50 uses the clean 50 gene list; the 22nd input is PAM50 plus clinical. Values are the mean and, in brackets, the range across the seven seeds. Ranked by mean V(BRACE) minus observed (`E_policy_ranking.csv`); the pooled column is from `E_seed_summary.csv`.

| Rank | Input | V(BRACE) minus observed | V(BRACE) minus treat none | No chemo arm (seeds defined) | Pooled RMST difference | Degenerate seeds |
|---|---|---|---|---|---|---|
| 1 | PAM50 plus clinical | 2.77 (−3.44 to 12.53) | 8.28 (3.95 to 15.08) | 32.50 (17.89 to 40.15), 7 | 18.77 (13.80 to 26.68) | 0 of 7 |
| 2 | Driver Foundation | 0.33 (−6.65 to 8.92) | 5.83 (1.14 to 11.60) | 5.33 (−2.56 to 17.20), 7 | 16.50 (14.72 to 18.76) | 0 |
| 3 | Reference: observed care | 0 | 5.50 (2.55 to 11.14) | | | |
| 4 | Reference: chemo for everyone | −1.78 (−9.69 to 8.15) | 3.72 (−6.36 to 10.70) | | | 7 |
| 5 | Driver Extended WGS | −2.24 (−8.41 to 11.68) | 3.26 (−2.73 to 14.23) | 24.85 (0.87 to 45.42), 7 | 18.82 (13.79 to 24.47) | 0 |
| 6 | Pearson correlation (100) | −2.61 (−9.29 to 8.22) | 2.89 (−5.96 to 10.77) | 8.61 (−13.76 to 16.12), 7 | 15.21 (10.36 to 18.18) | 0 |
| 7 | PAM50 (clean, 50) | −3.26 (−9.10 to 0.94) | 2.24 (−1.16 to 7.12) | 9.69 (3.13 to 13.73), 7 | 17.04 (13.33 to 22.35) | 0 |
| 8 | Driver Core | −3.45 (−12.95 to 7.18) | 2.05 (−4.33 to 9.73) | 7.98 (2.67 to 14.92), 7 | 16.60 (13.15 to 20.01) | 0 |
| 9 | Prolif Core CC | −4.16 (−7.41 to −1.11) | 1.34 (−1.23 to 5.04) | 14.48 (−16.05 to 40.50), 5 | 17.84 (11.88 to 25.66) | 0 |
| 10 | Control theory Type I, K out | −4.77 (−9.84 to −2.16) | 0.73 (−1.04 to 4.02) | 13.73 (−4.99 to 24.24), 6 | 18.00 (13.57 to 22.50) | 0 |
| 11 | TCGA BRCA Cohorts | −4.94 (−10.74 to −1.73) | 0.56 (−3.18 to 3.73) | 6.91 (2.06 to 17.06), 7 | 15.13 (11.99 to 21.81) | 0 |
| 12 | Blueprint | −4.99 (−11.00 to −2.51) | 0.52 (−0.81 to 3.67) | 0.07 (−8.77 to 11.37), 3 | 18.92 (13.77 to 26.09) | 2 |
| 13 | Interval Screen | −5.13 (−9.54 to −0.32) | 0.38 (−5.56 to 4.41) | 11.01 (−1.00 to 44.05), 7 | 17.03 (13.03 to 21.31) | 0 |
| 14 | TCGA BRCA All | −5.25 (−11.80 to −1.70) | 0.25 (−2.11 to 2.91) | 18.06 (8.53 to 39.59), 7 | 17.03 (13.98 to 19.19) | 0 |
| 15 | TCGA BRCA Subgroups | −5.47 (−11.14 to −2.34) | 0.03 (0.00 to 0.22) | none defined | 18.65 (13.77 to 24.35) | 6 |
| 16 | Reference: no chemo for anyone | −5.50 (−11.14 to −2.55) | 0 | | 18.66 (13.77 to 24.35) | 7 |
| 17 | EndoPredict | −5.50 (−11.14 to −2.55) | 0 | none defined | 18.66 (13.77 to 24.35) | 7 |
| 18 | Control theory Type II, K in | −5.51 (−11.70 to 1.40) | −0.01 (−2.84 to 3.95) | 10.81 (−15.48 to 37.64), 6 | 17.87 (12.35 to 23.33) | 0 |
| 19 | Control theory Type II, K out | −5.57 (−10.67 to −2.12) | −0.07 (−1.49 to 1.10) | 17.13 (−4.53 to 48.50), 7 | 16.90 (11.07 to 20.33) | 0 |
| 20 | Control theory Type I, K in | −5.73 (−11.47 to −2.95) | −0.22 (−0.91 to 0.72) | 19.71 (−11.98 to 49.87), 5 | 18.83 (13.76 to 25.37) | 0 |
| 21 | Control theory Type II, K | −5.73 (−11.42 to −2.41) | −0.23 (−0.92 to 0.92) | 8.22 (−3.46 to 20.17), 3 | 17.52 (13.18 to 21.29) | 0 |
| 22 | Index BCI | −5.96 (−11.05 to −3.05) | −0.46 (−2.57 to 0.30) | 17.24 (−10.75 to 63.74), 4 | 18.39 (14.08 to 22.52) | 1 |
| 23 | Oncotype | −6.27 (−10.47 to −1.75) | −0.77 (−7.15 to 4.79) | 2.51 (−15.83 to 27.40), 5 | 16.60 (9.61 to 22.89) | 1 |
| 24 | Control theory Type I, K | −6.28 (−11.35 to −2.55) | −0.78 (−4.06 to 0.16) | 6.17 (−16.01 to 19.87), 4 | 16.70 (12.93 to 21.07) | 1 |
| 25 | MammaPrint | −6.54 (−10.86 to −3.03) | −1.04 (−4.68 to 1.53) | 5.23 (−24.33 to 26.26), 6 | 18.35 (13.00 to 24.65) | 0 |

At seed 123 all 20 inputs other than PAM50 reproduce the pooled values of Table 3 to two decimals, and the arm counts for seed 150 are 122 chemo and 472 no chemo (Table S2 prints 222). The pooled comparison cannot rank inputs: recommending chemotherapy to nobody scores 18.66 months, which ties the top of that ranking, and the three identical entries of 17.12 in Table 3 (TCGA BRCA Subgroups, EndoPredict, Control theory Type I K) are exactly that degenerate policy at seed 123. On policy value, no input has a V(BRACE) minus observed range above zero across seeds, none of the 154 runs has a bootstrap interval above zero, only PAM50 plus clinical and Driver Foundation have a positive mean, and PAM50 alone ranks seventh with a mean of −3.26 months; 4 of 154 runs have an interval above zero for V(BRACE) minus treat none, about what chance alone would give.

Checklist items closed: 7 (matched seed comparison with uncertainty) and 33 (seed 150 arm counts).

---

## 7. Run G: recommendation intervals and abstention

Design (agreed in session): the seed in this pipeline only changes the split, so the seven E runs do not share test patients. G therefore keeps the seed 123 test set (594 patients) and refits CausalPFN (PAM50 plus clinical) on seven random 80% subsamples of the 1,386 training patients, one per S2 seed. A patient gets a confident recommendation only when all seven estimates have the same sign.

| Measure | Result |
|---|---|
| Confident recommendations | 378 of 594 (63.6%): 32 chemo, 346 no chemo; 216 abstain |
| Confident among C3 chemo recommendations / no chemo recommendations | 30.5% / 70.8% |
| Median per patient range of est_effect across the seven fits | 20.9 months |
| Agreement of the confident recommendations with C3 | 100% |
| No chemo arm (confident patients) | 23.37 (−0.79 to 47.53), p = 0.058; 286 Followed, 13 Not Followed |
| Chemo arm (confident patients) | −4.19 (−33.00 to 24.63), p = 0.78 |
| Pooled (confident patients) | 16.75 (4.02 to 29.48), p = 0.010 |
| Overlap policy value (104 confident patients in the overlap population) | V(BRACE) 115.6; minus observed −0.18 (−20.48 to 15.32); minus treat none 7.22 (−1.48 to 14.55) |

About two thirds of patients receive a stable recommendation, but most chemotherapy recommendations are unstable: only 32 of the 105 patients recommended chemotherapy in C3 keep that recommendation in all seven refits. Restricting to confident patients does not produce a policy value that beats observed care, and the no chemo arm comparison rests on 13 Not Followed patients. An abstention rule is easy to report and makes the uncertainty visible, but it does not by itself turn the recommendations into a demonstrated benefit.

Checklist items closed: 17 (recommendation intervals).

---

## 8. Run F: semi synthetic benchmark

Completed 11 October 2026, 23:26: 36 settings by 10 replicates (360 replicates), no failures, about 12 seconds per replicate. Script `F_semisynthetic.py`, per replicate results in `F_semisynthetic/reps/`, averages and standard deviations per setting and estimator in `F_semisynthetic/F_summary.csv`, log in `logs/F/F_run.log`.

Design: real METABRIC covariates (PAM50 clean plus clinical, seed 123 split, 1,386 training and 594 test rows) with simulated treatment and survival. Propensity e(x) = expit(a0 + z), where z is a standardised linear score of age, NPI, ER, ERBB2 and FOXA1 and a0 is set to give 50%, 20% or 10% treated (the slope is fixed at 1, so the degree of confounding is the same in every setting). Survival times follow a Weibull model in which chemotherapy changes log survival time by an amount depending on ESR1, MKI67 and age; 44% of test patients truly benefit. Censoring is uniform at 0%, 20%, 40% or 60% (achieved exactly), and 0, 50 or 200 pure noise features are added. Estimators: CausalPFN with the paper outcome min(U, 360) ("PFN paper"), CausalPFN with the IPCW outcome at 150 months ("PFN IPCW") and the grf causal survival forest at horizon 150 ("CSF"). Regret is the true value of the optimal policy minus the true value of the recommended policy, in months of RMST at 150; the optimal policy gains 10.6 months over the better of treating nobody or everybody.

| Setting | PEHE: PFN paper / PFN IPCW / CSF | Sign agreement: PFN paper / PFN IPCW / CSF | Regret: PFN paper / PFN IPCW / CSF |
|---|---|---|---|
| All 360 replicates | 31.4 / 21.2 / 22.5 | 0.78 / 0.74 / 0.71 | 3.0 / 3.9 / 4.6 |
| 50% treated | 28.8 / 19.5 / 21.5 | 0.84 / 0.82 / 0.80 | 1.8 / 2.4 / 1.8 |
| 20% treated | 31.9 / 20.4 / 21.6 | 0.80 / 0.75 / 0.73 | 2.4 / 3.5 / 4.0 |
| 10% treated | 33.6 / 23.7 / 24.6 | 0.70 / 0.66 / 0.59 | 4.9 / 5.8 / 7.9 |
| 0% censored | 44.1 / 17.9 / 22.7 | 0.83 / 0.81 / 0.70 | 1.6 / 2.2 / 4.7 |
| 20% censored | 33.9 / 19.7 / 22.4 | 0.81 / 0.78 / 0.72 | 2.4 / 2.9 / 4.3 |
| 40% censored | 24.7 / 21.7 / 22.3 | 0.77 / 0.74 / 0.71 | 3.1 / 4.0 / 4.5 |
| 60% censored | 23.0 / 25.3 / 22.8 | 0.71 / 0.65 / 0.70 | 4.9 / 6.6 / 4.9 |
| 0 noise features | 32.6 / 18.7 / 21.0 | 0.82 / 0.79 / 0.74 | 1.9 / 2.6 / 3.7 |
| 50 noise features | 31.4 / 22.2 / 22.6 | 0.75 / 0.72 / 0.71 | 3.7 / 4.5 / 4.6 |
| 200 noise features | 30.3 / 22.6 / 24.0 | 0.76 / 0.72 / 0.67 | 3.4 / 4.6 / 5.5 |
| Hardest corner (10% treated, 60% censored) | 25.1 / 27.4 / 24.8 | 0.61 / 0.57 / 0.58 | 7.1 / 8.6 / 8.2 |

Rows other than the first and last average over the other two factors. PEHE for "PFN paper" is not comparable with the other two, because that outcome targets a 360 month horizon while the truth is the effect on RMST at 150 months; its sign agreement and regret are comparable. The recommended policy beat the better constant policy in 96% (PFN paper), 94% (PFN IPCW) and 82% (CSF) of replicates.

When the data are simulated from a known model, CausalPFN does recover who benefits: with the paper outcome it has the highest sign agreement and the lowest regret overall, and it gives the lowest regret in 26 of the 36 settings (CSF in 9, PFN IPCW in 1), while PFN IPCW has lower PEHE than CSF in 24 of the 36 settings. All three estimators degrade as the treated share falls, censoring rises and noise features are added, and at 10% treated with 60% censoring, which is closest to METABRIC (about 23% treated), all three get only 57% to 61% of signs right and lose 7 to 9 of the 10.6 months available. This benchmark shows the estimator can work when its assumptions hold, but censoring here is independent of treatment and covariates and the effect is smooth, so it does not explain why the censoring aware outcomes in run D and the policy values in run E show no benefit on the real data.

Checklist items closed: 6 (semi synthetic benchmark).

---

## 9. Checklist items covered

| Item | Run | Status |
|---|---|---|
| 1 Clinical covariates in X | C | closed |
| 2 Censoring aware outcomes | D | closed |
| 6 Semi synthetic benchmark | F | closed |
| 7 Matched seed comparison | E | closed |
| 17 Recommendation intervals | G | closed |
| 24 Missing run_causalpfn.py | Setup | closed except the GitHub commit |
| 33 Seed 150 arm counts | E | closed (122 / 472) |
| 34 Checkpoint and UBE2T | Setup, B | closed |

## 10. Files

Scripts: `run_causalpfn.py`, `R/` (edited copies of the pipeline plus `revision_helpers.R`, `evaluate_run.R`, `csf_fit.R`), `build_clinical_covariates.py`, `censoring_outcomes.py`, `policy_value_for_run.py`, `append_summary.py`, `compare_to_saved.py`, `run_BCD.sh`, `run_EG.py`, `summarise_E.py`, `rank_E_policy.py`, `analyse_G.py`, `F_semisynthetic.py`. Results: `summary.csv`, `F_semisynthetic/F_summary.csv`, `E_seed_summary.csv`, `E_seed_long.csv`, `E_pam50_paired.csv`, `E_policy_ranking.csv`, `G_patient_spread.csv`, `G_summary.json`, per run CSV and JSON files, and `logs/`.
