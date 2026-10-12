# Preprocess evaluation dataset GSE20685

library(GEOquery)
library(dplyr)
library(stringr)

library(hgu133plus2.db)
library(AnnotationDbi)

options(GEOquery.inmemory.gpl = FALSE)
options(timeout = 600)

gse20685 <- getGEO("GSE20685", GSEMatrix = TRUE)[[1]]
clin <- pData(gse20685)
gene <- exprs(gse20685)

# rename columns we want to keep
clin <- clin %>%
  transmute(
    patientID = rownames(clin)
    , time = as.numeric(str_extract(`characteristics_ch1.4`, "[0-9.]+")) * 12 # transform to months
    , event = as.numeric(str_extract(`characteristics_ch1.3`, "[0-9]+"))
    , treatment = `adjuvant_chemotherapy:ch1`
  )

# histogram
hist(clin$time,
     breaks = 30,
     col = "skyblue",
     border = "white",
     main = "Distribution of Survival Time",
     xlab = "Time",
     ylab = "Frequency"
     )

# assume "unknown" treatment is "no"
table(clin$treatment)
clin <- clin %>%
  mutate(
    treatment = ifelse(treatment %in% c("yes"), 1, 0)
  )

# good, now it is binary
table(clin$treatment)

# format gene expressions
gene <- t(gene)
gene <- as.data.frame(gene)
gene$patientID <- rownames(gene)

# merge clinical and gene data
df <- clin %>%
  inner_join(gene, by = "patientID")

# good, all populated
colSums(is.na(df[, 1:4]))

# find annotation used for probe IDs
annotation(gse20685)

# map probe IDs to gene expressions
probe_ids <- colnames(df)[-(1:4)]

gene_symbols <- mapIds(
  hgu133plus2.db
  , keys = probe_ids
  , column = "SYMBOL"
  , keytype = "PROBEID"
  , multiVals = "first"
)

# rename columns
colnames(df)[-(1:4)] <- gene_symbols

# 9977 probes did not match to a gene
sum(is.na(gene_symbols))

# remove these
df <- df[, !is.na(colnames(df))]

# squish duplicate gene columns
  # (sometimes a probe contributes partly to multiple genes - aggregate and squish)
gene_data <- df[, -(1:4)]

gene_data_collapsed <- as.data.frame(
  t(rowsum(t(gene_data), group = colnames(gene_data)))
)

df <- cbind(df[, 1:4], gene_data_collapsed)

# save result
path = "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/GSE20685/"
write.csv(df, file=paste0(path, "dataset.csv"), row.names=FALSE)

