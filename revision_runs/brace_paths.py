"""Shared paths for the revision scripts (Python side; R/paths.R is the R mirror and builds the same layout).

BRACE_ROOT (environment variable) is the only path setting. Default: the parent folder of the repository.
  BRACE_ROOT/Processed_Datasets/   dataset_genes.csv and dataset_genes_and_clinical.csv
  BRACE_ROOT/Gene_Sets/            gene set files (e.g. 7_PAM50.csv)
  BRACE_ROOT/BRACE/Datasets/       analysis_extracts/ and metabric_evaluation_data_20260617.csv
Code and outputs live in the repository itself: revision_runs/ (this folder) and policy_value_analysis/.
BRACE_PYTHON (optional) is the Python interpreter used for child processes; default: the one running this script.
"""
import os
import sys
from pathlib import Path

REV = Path(__file__).resolve().parent              # revision_runs/ (scripts and per run outputs)
REPO = REV.parent
ROOT = Path(os.environ.get("BRACE_ROOT", str(REPO.parent))).expanduser().resolve()
DATA_DIR = ROOT / "Processed_Datasets"
GENE_SET_DIR = ROOT / "Gene_Sets"
DATASETS_DIR = ROOT / "BRACE" / "Datasets"
PV_SCRIPT = REPO / "policy_value_analysis" / "policy_value.py"
PYTHON = os.environ.get("BRACE_PYTHON", sys.executable)
