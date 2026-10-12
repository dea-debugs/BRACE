"""Build revision_runs/clinical_covariates.csv: one row per METABRIC patient, 10 raw clinical covariates.

From Processed_Datasets/dataset_genes_and_clinical.csv (already encoded by 1_create_nice_dataset.R):
  AGE_AT_DIAGNOSIS, LYMPH_NODES_EXAMINED_POSITIVE, NPI, ER_IHC (1 = positive), INFERRED_MENOPAUSAL_STATE (1 = post)
From cBioPortal datahub brca_metabric data_clinical_sample.txt, joined on PATIENT_ID:
  GRADE (1-3), TUMOR_SIZE (mm), TUMOR_STAGE (0-4), HER2 (1 = Positive), PR (1 = Positive)
Missing values are left as NA; median imputation with missing indicators is fitted on the training
partition inside master.R (clinical_features.R), never here.
"""

import pandas as pd

from brace_paths import DATA_DIR, REV

OUT = REV / "clinical_covariates.csv"

FROM_DATASET = ["AGE_AT_DIAGNOSIS", "LYMPH_NODES_EXAMINED_POSITIVE", "NPI", "ER_IHC", "INFERRED_MENOPAUSAL_STATE"]

ds = pd.read_csv(DATA_DIR / "dataset_genes_and_clinical.csv", usecols=["patientID"] + FROM_DATASET)
sam = pd.read_csv(REV / "cbioportal/data_clinical_sample.txt", sep="\t", comment="#")

sam = pd.DataFrame({
    "patientID": sam.PATIENT_ID,
    "GRADE": pd.to_numeric(sam.GRADE, errors="coerce"),
    "TUMOR_SIZE": pd.to_numeric(sam.TUMOR_SIZE, errors="coerce"),
    "TUMOR_STAGE": pd.to_numeric(sam.TUMOR_STAGE, errors="coerce"),
    "HER2": sam.HER2_STATUS.map({"Positive": 1.0, "Negative": 0.0}),
    "PR": sam.PR_STATUS.map({"Positive": 1.0, "Negative": 0.0}),
})
assert not sam.patientID.duplicated().any()

clin = ds.merge(sam, on="patientID", how="left", validate="one_to_one")
assert len(clin) == len(ds) == 1980, len(clin)
clin.to_csv(OUT, index=False)

print(f"wrote {OUT} ({len(clin)} patients)")
print("missing per column:", clin.isna().sum().to_dict())
