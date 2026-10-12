# METABRIC BRCA data
# Source: email from Sindy, 12/10/2025

memory.limit()
sessionInfo()

library(dplyr)
library(caret)
library(fastDummies)

# INPUT DATA ###################################################################
# define input folder path
in_path = "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/brca_metabric/"

# read in clinical data
clinical <- read.table(paste0(in_path, "data_clinical_patient.txt"), header = TRUE, sep = "\t")

# rename patientID column
names(clinical)[names(clinical) == "PATIENT_ID"] <- "patientID"

# read in gene expressions data
genes <- read.table(paste0(in_path, "data_mrna_illumina_microarray.txt"), header = TRUE, sep = "\t")

# transpose dataframe
genes <- as.data.frame(t(genes))

# make first row into header
header <- as.character(genes[1, ])
genes <- genes[-1, ]
genes <- genes[-1, ]
colnames(genes) <- header

# create patientID column
genes$patientID <- sapply(strsplit(rownames(genes), "-"), function(x) x[1])

# replace . with - for consistency with clinical dataframe
# use gsub to replace all occurrences
# . is a special character in regex, so use \\ to escape it
genes$patientID <- gsub("\\.", "-", genes$patientID)

# create master table
master_df <- merge(clinical, genes, by="patientID")

clinical_cols <- colnames(clinical)

# remove previous dataframes to save memory
rm(clinical)
rm(genes)
gc()

print_dims <- function(df) {
  dims <- dim(df)
  message("Rows: ", dims[1], ", Columns: ", dims[2])
}


# FEATURE ENGINEERING ##########################################################
# create survival columns
names(master_df)[names(master_df) == "OS_MONTHS"] <- "time"

master_df$event <- ifelse(master_df[["OS_STATUS"]] == "1:DECEASED", 1,
                          ifelse(master_df[["OS_STATUS"]] == "0:LIVING", 0, NA))

# Create new outcome column restricted survival time (RST) at 5 years [EDIT: this is 150 years]
tau <- 5 * 365

master_df$rst_5 <- mapply(
  function(time, event) {
    if (is.na(time) || is.na(event)) {
      return(NA)
    } else if (event == 1 && time <= tau) {
      return(time)
    } else {
      return(min(time, tau))
    }
  },
  master_df$time, master_df$event
)

master_df <- master_df %>% filter(!is.na(rst_5))

# select treatment column chemotherapy

# possible treatments
# CHEMOTHERAPY: YES | NO
# HORMONE_THERAPY: YES | NO
# RADIO_THERAPY: YES | NO
# BREAST_SURGERY: MASTECTOMY | BREAST CONSERVING

treatments <- c("CHEMOTHERAPY", "RADIO_THERAPY", "HORMONE_THERAPY", "BREAST_SURGERY")

treatment <- "CHEMOTHERAPY"

print_dims(master_df)

master_df$treatment <- ifelse(master_df[[treatment]] == "YES", 1, 0)
table(master_df$treatment)

master_df$all_no <- master_df$CHEMOTHERAPY == "NO" & master_df$RADIO_THERAPY == "NO" & master_df$HORMONE_THERAPY == "NO"
df_true <- master_df[master_df$all_no, ]

# remove specified list of columns
cols_to_remove <- c("OS_STATUS", "all_no")

df_clin <- master_df %>% select(-all_of(cols_to_remove))

clinical_features <- intersect(colnames(df_clin), clinical_cols)
categorical_cols <- clinical_features[sapply(df_clin[clinical_features], \(x) is.factor(x) || is.character(x))]
cat("Categorical clinical columns to be one-hot encoded:\n")
print(categorical_cols)

cat("\nPreviewing first 10 unique values for each categorical column:\n")
for (col in categorical_cols) {
  cat("\n----", col, "----\n")
  print(unique(df_clin[[col]])[1:10])
}

df_encoded <- df_clin

# binarise
to_binary <- c("ER_IHC"
               , "INFERRED_MENOPAUSAL_STATE"
               , "LATERALITY"
               , "BREAST_SURGERY"
               , "RFS_STATUS"
               )

# clean up
for (col in to_binary) {
  df_encoded[[col]] <- trimws(df_encoded[[col]])
  df_encoded[[col]] <- ifelse(is.na(df_encoded[[col]]) | df_encoded[[col]] == "", NA, df_encoded[[col]])
}

# custom binary mappings
df_encoded$ER_IHC <- ifelse(df_encoded$ER_IHC == "Positve", 1,
                            ifelse(df_encoded$ER_IHC == "Negative", 0, NA))

df_encoded$INFERRED_MENOPAUSAL_STATE <- ifelse(df_encoded$INFERRED_MENOPAUSAL_STATE == "Post", 1,
                                               ifelse(df_encoded$INFERRED_MENOPAUSAL_STATE == "Pre", 0, NA))

df_encoded$LATERALITY <- ifelse(df_encoded$LATERALITY == "Right", 1,
                                ifelse(df_encoded$LATERALITY == "Left", 0, NA))

df_encoded$BREAST_SURGERY <- ifelse(df_encoded$BREAST_SURGERY == "MASTECTOMY", 1,
                                    ifelse(df_encoded$BREAST_SURGERY == "BREAST CONSERVING", 0, NA))

df_encoded$RFS_STATUS <- ifelse(df_encoded$RFS_STATUS == "1:Recurred", 1,
                                ifelse(df_encoded$RFS_STATUS == "0:Not Recurred", 0, NA))

# one-hot encode
to_encode <- c("CELLULARITY"
               , "HER2_SNP6"
               , "INTCLUST"
               , "THREEGENE"
               , "HISTOLOGICAL_SUBTYPE"
               , "CLAUDIN_SUBTYPE"
               )

df_encoded <- fastDummies::dummy_cols(
  df_encoded,
  select_columns = to_encode,
  remove_selected_columns = TRUE,
  ignore_na = TRUE
)

df_encoded <- df_encoded %>% select(-c(VITAL_STATUS, CHEMOTHERAPY, HORMONE_THERAPY, SEX, RADIO_THERAPY, RFS_MONTHS, RFS_STATUS, COHORT))

write.csv(df_encoded, file=paste0(in_path, "dataset_genes_and_clinical.csv"), row.names=FALSE)

# remove previous dataframes to save memory
rm(df_clin)
rm(df_encoded)
gc()

clinical_cols <- setdiff(clinical_cols, c("OS_MONTHS", "patientID"))
cols_to_remove <- unique(c(cols_to_remove, treatment, clinical_cols))

master_df <- master_df %>% select(-all_of(cols_to_remove))

# plot distributions
png(paste0(in_path, "results/", "chemo_histogram.png"), width = 800, height = 600)
hist(master_df$treatment, 
     main = "Histogram of treatment", 
     xlab = "Treatment", 
     col = "skyblue", 
     border = "white")
dev.off()

png(paste0(in_path, "results/", "rst_histogram.png"), width = 800, height = 600)
hist(master_df$rst_5, 
     main = "Histogram of 5-yr RST", 
     xlab = "Outcome", 
     col = "skyblue", 
     border = "white")
dev.off()

# remove no-longer-necessary variables from memory
rm(header, treatments, tau, cols_to_remove, treatment)
gc()

# save dataset to save time next time
write.csv(master_df, file=paste0(in_path, "dataset_genes.csv"), row.names=FALSE)
