# Shared paths for the revision R scripts (R side; ../brace_paths.py is the Python mirror and builds the same layout).
# BRACE_ROOT (environment variable) is the only path setting. Default: the parent folder of the repository.
#   BRACE_ROOT/Processed_Datasets/   dataset_genes.csv and dataset_genes_and_clinical.csv
#   BRACE_ROOT/Gene_Sets/            gene set files (e.g. 7_PAM50.csv)
# Code and outputs live in the repository itself: revision_runs/ (REVISION_DIR) and revision_runs/R/ (this folder).
# BRACE_PYTHON (optional) is the Python interpreter of the brace conda environment; default: python on the PATH.

.brace_r_dir <- function() {
  for (i in rev(seq_len(sys.nframe()))) {  # innermost source() call that has a file name, i.e. this file
    f <- sys.frame(i)$ofile
    if (!is.null(f) && basename(f) == "paths.R") return (dirname(normalizePath(f)))
  }
  normalizePath(getwd())                   # fallback: scripts are run from revision_runs/R
}

R_DIR        <- .brace_r_dir()
REVISION_DIR <- paste0(dirname(R_DIR), "/")
BRACE_ROOT   <- normalizePath(Sys.getenv("BRACE_ROOT", unset = dirname(dirname(dirname(R_DIR)))),  # R/ -> revision_runs/ -> repo -> parent
                              mustWork = FALSE)
DATA_DIR     <- paste0(BRACE_ROOT, "/Processed_Datasets/")
GENE_SET_DIR <- paste0(BRACE_ROOT, "/Gene_Sets/")
CLINICAL_FILE <- paste0(REVISION_DIR, "clinical_covariates.csv")

CAUSALPFN_PYTHON <- Sys.getenv("BRACE_PYTHON", unset = unname(Sys.which("python")))
CAUSALPFN_SCRIPT <- paste0(REVISION_DIR, "run_causalpfn.py")
CAUSALPFN_IO_DIR <- paste0(REVISION_DIR, "causalpfn_io")
CENSORING_SCRIPT <- paste0(REVISION_DIR, "censoring_outcomes.py")
