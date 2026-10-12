# evaluation_20260617

library(dplyr)

# set up file paths
in_path <- "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/brca_metabric/"
out_dir <- "C:/Users/parki/OneDrive/Desktop/Honours/HTE_results/"
dataset <- "gse_meta_evaluation"

# read in datasets
predictions <- read.csv(
  paste0(out_dir, dataset, "_20260617.csv")
)

clinical <- read.csv(
  paste0(in_path, "dataset_genes_and_clinical.csv")
)

# join datasets
df <- predictions %>%
  left_join(clinical, by = "patientID")

# format subtypes
subtype_cols <- grep("^CLAUDIN_SUBTYPE_", names(df), value = TRUE)

df <- df %>%
  mutate(
    subtype = sub("^CLAUDIN_SUBTYPE_", "",
                  subtype_cols[max.col(across(all_of(subtype_cols)),
                                       ties.method = "first")])
  )

df <- df %>%
  select(-all_of(subtype_cols))

# summary stats
n_patients <- nrow(df)
cat("Total number of patients:", n_patients, "\n\n")

df %>%
  count(subtype) %>%
  mutate(
    percentage = round(100 * n / sum(n), 2)
  ) %>%
  arrange(desc(n))

# output dataset
write.csv(
  df,
  file = paste0(out_dir, dataset, "_data_20260617", ".csv"),
  row.names = FALSE
)
