# -----------------------------------------------
# Inputs: X_train, X_test.
# Methods are to run on X_train only, to avoid data leakage.
# Outputs: X_train, X_test - including selected features only.
# -----------------------------------------------

library(FSelectorRcpp)
source("utils.R")
if (!exists("GENE_SET_DIR")) source("paths.R")

fs_mutual_information <- function(X_train, X_test, y_train, n = 100) {
  
  # format data as required
  df <- data.frame(X_train, y = y_train$rst_30)
  
  # calculate information
  mi <- information_gain(y ~ ., df)
  
  # rank by importance
  mi <- mi[order(mi$importance, decreasing = TRUE), ]
  
  # select top features
  selected <- mi$attributes[1:min(n, nrow(mi))]
  
  return (list(X_train = X_train[, selected, drop = FALSE]
               , X_test  = X_test[,  selected, drop = FALSE]))
}

fs_pearson_correlation <- function(X_train, X_test, y_train, n = 100) {
  
  # calculate correlation
  cors <- sapply(X_train, function(x) {
    cor(x, y_train$rst_30, use = "complete.obs")
  })
  
  # rank by absolute correlation
  ranked <- sort(abs(cors), decreasing = TRUE)
  
  # select top features
  selected <- names(ranked)[1:min(n, length(ranked))]
  
  return (list(X_train = X_train[, selected, drop = FALSE]
          , X_test  = X_test[,  selected, drop = FALSE]))
}

fs_gene_set <- function(X_train, X_test, set, n = 100,
                        path = GENE_SET_DIR,
                        gene_file = NULL) {
  # Options for set:
  # 1, 2, 3, 4A, 4B, 4C, 5, 6, 7, 8, 9, 10, 11, 12
  # REVISION EDIT: local gene set folder; gene_file reads an explicit list instead (e.g. 7_PAM50_clean.csv)

  if (is.null(gene_file)) {
    features <- utils_find_file_with_prefix(path=path, prefix=paste0(set, "_"))
  } else {
    features <- read.csv(gene_file)
  }

  # select features
  selected <- features$Gene
  
  missing <- setdiff(selected, colnames(X_train))
  print(paste0("Missing: ", missing))
        
  selected <- intersect(selected, colnames(X_train))
        
  
  return (list(X_train = X_train[, selected, drop = FALSE]
               , X_test  = X_test[,  selected, drop = FALSE]))
  
}

fs_control_theory <- function(X_train, X_test, type, direction, n = 100,
                              path = GENE_SET_DIR) {
# Options for type:
  # I, II
  # Options for direction:
  # None, in, out
  
  prefix <- paste0("filtered_ppi_network_Type", type, "_ranked_by_k")
  
  add_suffix <- function(x) {
    if (is.null(x)) return("")
    paste0("_", x)
  }
  
  prefix <- paste0(prefix, add_suffix(direction))

  # REVISION EDIT: local gene set folder (path argument)
features <- utils_find_file_with_prefix(path=path, prefix=prefix)
  
  # order features by ranking
  features <- features[order(features$rank), ]
  
  # select top features
  selected <- features$Gene[1:n]
  
  missing <- setdiff(selected, colnames(X_train))
  print(paste0("Missing: ", missing))
  
  selected <- intersect(selected, colnames(X_train))
  
  return (list(X_train = X_train[, selected, drop = FALSE]
               , X_test  = X_test[,  selected, drop = FALSE]))
}

fs_combo <- function(X_train, X_test, set) {
  # Options for set:
  # ct_pam50, ct_pam50_100
  
  path <- "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/brca_metabric/gene_sets/"
  features <- read.csv(paste0(path, set, '.csv'))
  
  # select features
  selected <- features$Gene
  
  missing <- setdiff(selected, colnames(X_train))
  print(paste0("Missing: ", missing))
  
  selected <- intersect(selected, colnames(X_train))
  
  return (list(X_train = X_train[, selected, drop = FALSE]
               , X_test  = X_test[,  selected, drop = FALSE]))
  
}