# -----------------------------------------------
# Goal 1: To select a Heterogeneous Treatment Effect (HTE) method.
# Goal 2: To select a feature selection (FS) method.

# Definitions of methods ------------------------
# Log-rank test
# Null hypothesis:  All groups are the same
# Alternate hypothesis:  At least one group differs
# If p-value is small, that's good for our purposes.
# This test is to check the quality of the partitions, it is not a causal estimand itself.

# Interpretation and consistency of methods -----
# Causal tree outputs explicit partitions, via leaves.
# Other methods (CF, causalPFN) output continuous results and do not naturally define explicit groups.

# Concordance index:  Directly comparable across models.

# P-value:  Need to construct groups to enable comparison across models.
# Could implement quantile-based grouping, like splitting into quartiles.
# But then, the P-value may not have much meaning unless the patients coincidentally cluster into four distinct separate groups.
# If they naturally grouped into different numbers, that would be obscured here.
# Choice of K is arbitrary but mandatory.
# However, we're not testing for specific clusters.  We're testing for risk stratification.
# i.e., Does the model induce an ordering of patients such that survival differs across the ordering?

# We'd probably want to do this for CTs as well, rather than using its own native leaves, for consistency.
# -----------------------------------------------

# New edit:  18/06/2026
# To run the balanced version, use:
#   - training "metabric"
#   - eval_type "internal"
#   - eval_set "metabric"
#   - change line 79 from "dataset_genes.csv" to "dataset_genes_balanced.csv"

# -- REVISION EDIT: run settings --
# A driver script may set any of these before source("master.R"); otherwise the defaults below apply.
# Defaults reproduce run A (PAM50 from the original list, seed 123, outcome rst_30).
if (!exists("RUN_LABEL"))  RUN_LABEL  <- "A_baseline_pam50_seed123"
if (!exists("SEED"))       SEED       <- 123
if (!exists("GENE_SET"))   GENE_SET   <- "7"      # prefix in GENE_SET_DIR, as in the original code
if (!exists("GENE_FILE"))  GENE_FILE  <- NULL     # optional explicit gene list file (overrides GENE_SET)
if (!exists("OUTCOME"))    OUTCOME    <- "rst_30"  # or "ipcw150", "pseudo_rmst150" (computed on train)
if (!exists("USE_CLINICAL")) USE_CLINICAL <- FALSE  # add the 10 clinical covariates (run C/D)
if (!exists("FS_CHOICE"))  FS_CHOICE  <- NULL     # "correlation" or "control_theory" (run E); NULL = gene set
if (!exists("CT_TYPE"))    CT_TYPE    <- "I"      # control theory: "I" or "II"
if (!exists("CT_DIRECTION")) CT_DIRECTION <- ""   # control theory: "" (ranked by K), "in", "out"
if (!exists("TRAIN_SUBSAMPLE")) TRAIN_SUBSAMPLE <- NULL  # run G: fraction of the training partition kept
if (!exists("SUBSAMPLE_SEED"))  SUBSAMPLE_SEED  <- NULL  # run G: seed for that subsample (split still uses SEED)
CLINICAL_FILE <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/clinical_covariates.csv"
REVISION_DIR <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/"
DATA_DIR     <- "/Users/sindypinero/Downloads/4_Deanna/Processed_Datasets/"
GENE_SET_DIR <- "/Users/sindypinero/Downloads/4_Deanna/Gene_Sets/"
run_started  <- Sys.time()

# -- Import libraries --

library(dplyr)
library(survival)
library(causalTree)
library(grf)
library(randomForestSRC)
library(zeallot)
library(caret)

# -- Import custom functions --
source("utils.R")
source("feature_selection.R")
source("hte.R")
source("metrics.R")
source("revision_helpers.R")

set.seed(SEED)

out_dir <- REVISION_DIR
dir.create(out_dir, showWarnings = FALSE)
pdf(paste0(out_dir, RUN_LABEL, "_plots.pdf"))

# -- Data setup --

training <- "metabric"

# select subtype (if applicable)
# options: LumA, LumB, Her2, claudin.low, Basal, Normal, NC
subtype <- ""

# select evaluation dataset
# can be either internal or external
eval_type <- "internal"
eval_set <- "metabric"

# Get training and testing data
in_path <- DATA_DIR

if (training == "gse") {
  df <- read.csv(paste0("C:/Users/parki/OneDrive/Desktop/Honours/Datasets/GSE20685/", "dataset.csv"))

} else {
  
  if (subtype == "") {
    # REVISION EDIT: read revision_runs/cache/dataset_genes.rds when it matches the CSV (size + mtime)
    df <- read_dataset_cached(paste0(in_path, "dataset_genes.csv"), paste0(REVISION_DIR, "cache/"))

  } else {
    df <- read.csv(paste0(in_path, "dataset_genes_and_clinical.csv"))
    
    # filter for subtype
    colname <- grep(paste0("^CLAUDIN_SUBTYPE_.*", subtype, "$"), names(df), value=TRUE)
    df <- df %>% filter(.data[[colname]] == 1)
    
    # remove unwanted columns
    remove <- grep("^CLAUDIN_SUBTYPE_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    remove <- grep("^HER2_SNP6_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    remove <- grep("^INTCLUST_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    remove <- grep("^THREEGENE_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    remove <- grep("^HISTOLOGICAL_SUBTYPE_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    remove <- grep("^CELLULARITY_", names(df), value = TRUE)
    df <- df %>% dplyr::select(-all_of(remove))
    
    # drop unnecessary clinical columns
    to_drop <- c("LYMPH_NODES_EXAMINED_POSITIVE","NPI",
                "ER_IHC","INFERRED_MENOPAUSAL_STATE",
                "AGE_AT_DIAGNOSIS","LATERALITY"
                ,"BREAST_SURGERY"
              )
    
    df <- df %>% dplyr::select(-all_of(to_drop))
    
  }
}

# Create new outcome column restricted survival time (RST) at 30 years
tau <- 360

df$rst_30 <- mapply(
  function(time, event) {
    if (is.na(time) || is.na(event)) {
      return(NA)
    } else if (event == 1 && time <= tau) {
      return(time)
    } else {
      return(min(time, tau))
    }
  },
  df$time, df$event
)

tau <- 5

df$rst_5 <- mapply(
  function(time, event) {
    if (is.na(time) || is.na(event)) {
      return(NA)
    } else if (event == 1 && time <= tau) {
      return(time)
    } else {
      return(min(time, tau))
    }
  },
  df$time, df$event
)
# END EDIT

if (training != "gse") {
  c(df_train, df_test) %<-% utils_train_test_split(df)#, train_size = 0.3)

  # REVISION EDIT (run G): refit on a random subsample of the training partition, test set unchanged
  if (!is.null(TRAIN_SUBSAMPLE)) {
    set.seed(SUBSAMPLE_SEED)
    df_train <- df_train[sort(sample(nrow(df_train), round(TRAIN_SUBSAMPLE * nrow(df_train)))), ]
  }

  # Remove from memory
  rm(df)
  gc()
  
  # Training features
  X_train <- df_train %>% dplyr::select(-c(patientID, treatment, time, event, rst_5, rst_30))
  y_train <- df_train[, c("treatment", "time", "event", "rst_5", "rst_30")]
  
  # Testing features
  X_test <- df_test %>% dplyr::select(-c(patientID, treatment, time, event, rst_5, rst_30))
  y_test <- df_test[, c("treatment", "time", "event", "rst_5", "rst_30")]
  
  # Patient info
  df_patient <- df_test[, c("patientID", "treatment", "time", "event", "rst_5", "rst_30")]
  
} else {
  df_train <- df
  X_train <- df %>% dplyr::select(-c(patientID, treatment, time, event, rst_5, rst_30))
  y_train <- df[, c("treatment", "time", "event", "rst_5", "rst_30")]
  
  # we won't actually use this X_test, we just need it to exist for feature selection
  X_test <- df %>% dplyr::select(-c(patientID, treatment, time, event, rst_5, rst_30))
}
# -----------------------------------------------
#
# Here, we perform feature selection.
#
# Inputs differ by function, so check feature_selection.R for more visibility.
#
# Options: mutual_information, correlation, control_theory (I, II | None, In, Out), gene_set (1, 2, 3, 4A, 4B, 4C, 5, 6, 7, 8, 9, 10, 11, 12)
# ct_pam50, ct_pam50_60, ct_pam50_100 - this is type I, ranked by K in
# ct_1_k_pam50, ct_1_k_pam50_100
# ct_1_kout_pam50, ct_1_kout_pam50_100
# ct_2_kin_pam50, ct_2_kin_pam50_100
# ct_2_k_pam50, ct_2_k_pam50_100
# ct_2_kout_pam50, ct_2_kout_pam50_100
#
#
# -----------------------------------------------

FS_METHOD <- if (!is.null(FS_CHOICE)) FS_CHOICE else if (identical(GENE_SET, "none")) "none" else "gene_set"
n <- 100

# REVISION EDIT: no genes (clinical only variant)
if (FS_METHOD == "none") {
  X_train <- X_train[, 0, drop = FALSE]
  X_test  <- X_test[, 0, drop = FALSE]
}

if (FS_METHOD == "mutual_information") {
  c(X_train, X_test) %<-% fs_mutual_information(X_train, X_test, y_train, n)
}

if (FS_METHOD == "correlation") {
  c(X_train, X_test) %<-% fs_pearson_correlation(X_train, X_test, y_train, n)
}

if (FS_METHOD == "gene_set") {
  set <- GENE_SET

  c(X_train, X_test) %<-% fs_gene_set(X_train, X_test, set, n, path = GENE_SET_DIR, gene_file = GENE_FILE)
}
  
if (FS_METHOD == "control_theory") {
  type <- CT_TYPE
  direction <- if (CT_DIRECTION == "") NULL else CT_DIRECTION

  c(X_train, X_test) %<-% fs_control_theory(X_train, X_test, type, direction, n, path = GENE_SET_DIR)
}

if (FS_METHOD == "combo") {
  set <- "ct_2_kout_pam50_100"
  
  c(X_train, X_test) %<-% fs_combo(X_train, X_test, set)
}

# REVISION EDIT: clinical covariates, imputation fitted on the training partition
if (USE_CLINICAL) {
  c(X_train, X_test) %<-% add_clinical_features(X_train, X_test, df_train$patientID, df_test$patientID, CLINICAL_FILE)
}

stopifnot(identical(colnames(X_train), colnames(X_test)))

# SWAP TO VALIDATION SET IF REQUESTED!!!
if (eval_type == "external") {
  df_test <- prepare_eval(eval_set, colnames(X_train))
  
  # old version
  #feature_cols <- colnames(X_train)
  #X_test <- df_test %>%
  #  dplyr::select(all_of(feature_cols))
  
  # new version
  feature_cols <- intersect(
    colnames(X_train),
    colnames(df_test)[-(1:4)]
  )
  
  X_train <- X_train %>% select(all_of(feature_cols))
  X_test  <- df_test %>% select(all_of(feature_cols))
  df_patient <- df_test[, c("patientID", "treatment", "time", "event", "rst_5", "rst_30")] # this is basically y_test
  
  # safety net check
  stopifnot(all(feature_cols %in% colnames(df_test)))
  stopifnot(identical(colnames(X_train), colnames(X_test)))
}

# -----------------------------------------------
#
# Here, we train the model.  The inputs are the X features, Treatment, and Outcome.
# The inputted outcome is either [Time and Event] or [RST].
# We then use the trained model to make predictions on the test set.
# The predictions are made only using the X features.
# Each of these can be rewritten as a function.

# The output of this section is estimated treatment effects for the test population.
#
# Options: random, causal_tree, causal_forest, causal_pfn
#
# -----------------------------------------------

# -- Define and run HTE method --

HTE_METHOD <- "causal_pfn"
risk_score <- NULL

if (HTE_METHOD == "random") {
  est_effect <- runif(nrow(X_test), min = -1, max = 1)
}

if (HTE_METHOD == "causal_pfn") {
  if (!OUTCOME %in% colnames(df_train)) df_train[[OUTCOME]] <- censoring_aware_outcome(df_train, OUTCOME)
  est_effect <- hte_causal_pfn(df_train, X_train, X_test, outcome = OUTCOME)
causalpfn_seconds <- attr(est_effect, "causalpfn_seconds")
  est_effect <- as.numeric(est_effect)
}

# Remove from memory
rm(df_train, df_test)
gc()

if (HTE_METHOD == "causal_tree") {
  est_effect <- hte_causal_tree(y_train, X_train, X_test)
}

if (HTE_METHOD == "causal_forest") {
  est_effect <- hte_causal_forest(y_train, X_train, X_test)
}

features_used <- colnames(X_train)

# Remove from memory
rm(X_test)
gc()

# Ensure risk_score size aligns with y_test size
if (length(est_effect) != nrow(df_patient)) {
  warning("Estimated Treatment Effect length does not match test data.")
}

# -----------------------------------------------
#
# Here, we make recommendations for each patient on whether they should take the treatment.
# We do this by setting a threshold.  This is a place we can apply tuning.
# We can set a decimal number, or we can use quartiles.
#
# -----------------------------------------------

# Safety net to check correct number of rows
stopifnot(nrow(df_patient) == length(est_effect))
stopifnot(!any(duplicated(df_patient$patientID)))

hist(est_effect,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Distribution of estimated treatment effect",
     xlab = "Values")

df_patient <- df_patient %>% mutate(est_effect = est_effect)

df_patient <- utils_recommend(df_patient, threshold = 0)

metric_coverage(df_patient$treatment, df_patient$recommendation)

#df_patient <- df_patient %>% filter(recommendation == 1)

# Analyse performance
metric_significance(df_patient)
metric_significance(df_patient, cutoff = 175)
metric_significance(df_patient, cutoff = 150)

#metric_hte_concordance(df_patient)
#metric_hte_concordance(df_patient, cutoff = 150)
#metric_hte_concordance(df_patient, cutoff = 175)

metric_survival_plot(df_patient)
metric_survival_plot(df_patient, lines = c(150, 175))
#metric_survival_plot(df_patient, lines = c(140))

library(survRM2)

df_patient$followed <- ifelse(df_patient$followed == "Y", 1,
                              ifelse(df_patient$followed == "N", 0, NA))

difference <- rmst2(
  time   = df_patient$time,
  status = df_patient$event,
  arm    = df_patient$followed,
  tau    = 175
)

# RMST results table
difference

print(paste0("Estimated survival difference: ", difference$unadjusted.result[1], " months"))
print(paste0("Survival difference P value: ", difference$unadjusted.result[10]))

# EXPORT FOR FURTHER EVALUATION
evaluation_export <- df_patient %>%
  dplyr::select(
    patientID,
    time,
    event,
    rst_30,
    est_effect,
    treatment,
    recommendation,
    followed
  )

# REVISION EDIT: per-run CSV and JSON in revision_runs/
run_export <- df_patient %>%
  mutate(run_label = RUN_LABEL) %>%
  dplyr::select(patientID, treatment, time, event, est_effect, recommendation, run_label)

write.csv(
  run_export,
  file = paste0(out_dir, RUN_LABEL, ".csv"),
  row.names = FALSE
)

jsonlite::write_json(
  list(
    run_label = RUN_LABEL,
    seed = SEED,
    training = training,
    eval_type = eval_type,
    eval_set = eval_set,
    data_file = paste0(in_path, "dataset_genes.csv"),
    split = "utils_train_test_split, train_size = 0.7, after set.seed(seed)",
    n_train = length(y_train$treatment),
    n_test = nrow(df_patient),
    fs_method = FS_METHOD,
    ct_type = if (FS_METHOD == "control_theory") CT_TYPE else NULL,
    ct_direction = if (FS_METHOD == "control_theory") CT_DIRECTION else NULL,
    train_subsample = TRAIN_SUBSAMPLE,
    subsample_seed = SUBSAMPLE_SEED,
gene_set = if (is.null(GENE_FILE)) GENE_SET else NULL,
    gene_file = GENE_FILE,
    n_features = length(features_used),
    features = features_used,
    hte_method = HTE_METHOD,
    outcome = OUTCOME,
    use_clinical = USE_CLINICAL,
    clinical_file = if (USE_CLINICAL) CLINICAL_FILE else NULL,
    clinical_imputation = if (USE_CLINICAL) "median + missing indicators, fitted on training partition" else NULL,
recommendation_threshold = 0,
    n_recommended_chemo = sum(df_patient$recommendation == 1),
    n_recommended_no_chemo = sum(df_patient$recommendation == 0),
    rmst_followed_vs_not_tau175 = unname(difference$unadjusted.result[1, 1]),
    rmst_followed_vs_not_tau175_p = unname(difference$unadjusted.result[1, 4]),
    causalpfn_seconds = causalpfn_seconds,
    causalpfn_python = CAUSALPFN_PYTHON,
    causalpfn_script = CAUSALPFN_SCRIPT,
    r_version = R.version.string,
    started = format(run_started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  ),
  paste0(out_dir, RUN_LABEL, ".json"),
  auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null"
)

dev.off()