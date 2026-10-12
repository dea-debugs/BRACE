# REVISION helpers sourced by master.R (revision_runs copy only).

# Clinical covariates for run C/D: join on patientID, then median imputation plus missing indicators
# fitted on the training partition only (indicator added for every column with missingness in train).
add_clinical_features <- function(X_train, X_test, ids_train, ids_test, clinical_file) {
  clin <- read.csv(clinical_file)
  clin_cols <- setdiff(colnames(clin), "patientID")

  c_train <- clin[match(ids_train, clin$patientID), clin_cols, drop = FALSE]
  c_test  <- clin[match(ids_test,  clin$patientID), clin_cols, drop = FALSE]
  stopifnot(!anyNA(match(ids_train, clin$patientID)), !anyNA(match(ids_test, clin$patientID)))

  for (col in clin_cols) {
    med <- median(c_train[[col]], na.rm = TRUE)
    if (anyNA(c_train[[col]])) {
      c_train[[paste0(col, "_miss")]] <- as.numeric(is.na(c_train[[col]]))
      c_test[[paste0(col, "_miss")]]  <- as.numeric(is.na(c_test[[col]]))
    }
    c_train[[col]][is.na(c_train[[col]])] <- med
    c_test[[col]][is.na(c_test[[col]])]   <- med
  }

  return (list(X_train = cbind(X_train, c_train), X_test = cbind(X_test, c_test)))
}

# Censoring aware outcomes for run D, computed on the training partition by censoring_outcomes.py.
# outcome: "ipcw150" (IPCW pseudo outcome) or "pseudo_rmst150" (jackknife pseudo observations of RMST).
censoring_aware_outcome <- function(df_train, outcome, io_dir = file.path(CAUSALPFN_IO_DIR, RUN_LABEL)) {
  dir.create(io_dir, showWarnings = FALSE, recursive = TRUE)
  method  <- sub("[0-9]+$", "", outcome)
  horizon <- as.numeric(sub("^[a-z_]+", "", outcome))
  in_file  <- file.path(io_dir, "outcome_input.csv")
  out_file <- file.path(io_dir, "outcome_output.csv")
  unlink(out_file)
  write.csv(df_train[, c("patientID", "time", "event", "treatment")], in_file, row.names = FALSE)

  script <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/censoring_outcomes.py"
  status <- system(paste(shQuote(CAUSALPFN_PYTHON), shQuote(script), shQuote(in_file), shQuote(out_file),
                         method, horizon))
  if (status != 0 || !file.exists(out_file)) stop("censoring_outcomes.py failed")

  y <- read.csv(out_file)
  stopifnot(identical(y$patientID, df_train$patientID))
  file.copy(sub("\\.csv$", ".json", out_file), paste0(REVISION_DIR, RUN_LABEL, "_outcome.json"), overwrite = TRUE)
  return (y$Y)
}

# Cached read of the large expression CSV: identical data frame to read.csv(csv_file), stored as an
# uncompressed RDS together with the CSV's size and mtime; rebuilt automatically if the CSV changes.
read_dataset_cached <- function(csv_file, cache_dir) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  base <- sub("\\.csv$", "", basename(csv_file))
  rds_file <- file.path(cache_dir, paste0(base, ".rds"))
  stamp_file <- file.path(cache_dir, paste0(base, ".source.txt"))
  info <- file.info(csv_file)
  stamp <- paste(normalizePath(csv_file), info$size, format(info$mtime, "%Y-%m-%d %H:%M:%OS3"))
  if (file.exists(rds_file) && file.exists(stamp_file) && identical(readLines(stamp_file), stamp)) {
    message("Reading cached ", rds_file)
    return (readRDS(rds_file))
  }
  message("Cache missing or stale, reading ", csv_file)
  df <- read.csv(csv_file)
  saveRDS(df, rds_file, compress = FALSE)
  writeLines(stamp, stamp_file)
  return (df)
}
