utils_train_test_split <- function(df, train_size = 0.7) {
  
  # Train/test split - default 70% training, 30% testing
  train_idx <- sample(seq_len(nrow(df)), size = train_size * nrow(df))
  df_train <- df[train_idx, ]
  df_test <- df[-train_idx, ]
  
  return (list(df_train, df_test))
}

utils_recommend <- function(df_patient, threshold = 0) {
  
  # Make recommendation for each patient
  df_patient$recommendation <- ifelse(df_patient$est_effect > threshold, 1, 0)
  
  # Did the patient follow our recommendation?
  df_patient$followed <- ifelse(df_patient$treatment == df_patient$recommendation, 'Y', 'N')
  
  return (df_patient)
}

utils_find_file_with_prefix <- function(path, prefix) {
  
  # locate correct file
  files <- list.files(
    path    = path,
    pattern = paste0("^", prefix),
    full.names = TRUE
  )
  
  print(path)
  print(prefix)
  print(files)
  
  # expect only one file
  stopifnot(length(files) >= 1)
  
  if (length(files) > 1) {
    files <- files[which.min(nchar(basename(files)))]
  }
  
  file <- files[1]
  ext  <- tolower(tools::file_ext(file))
  
  # read it in
  ext <- tools::file_ext(file)
  
  features <- switch(
    ext,
    csv = read.csv(file),
    txt = read.table(file, header = TRUE),
    tsv = read.table(file, header = TRUE, sep = "\t"),
    stop("Unsupported file type: ", ext)
  )
  
  return (features)
}

prepare_eval <- function(dataset, training_cols) {
  # input colnames(X_train)
  
  # read in evaluation dataset
  if (dataset == "gse") {
    path = "C:/Users/parki/OneDrive/Desktop/Honours/Datasets/GSE20685/"
    eval <- read.csv(paste0(path, "dataset.csv"))
  
  } else {
    eval <- read.csv(paste0("C:/Users/parki/OneDrive/Desktop/Honours/Datasets/brca_metabric/", "dataset_genes.csv"))
  
  }
  
  non_feature_cols <- c("patientID", "time", "event", "treatment", "rst_5", "rst_30")
  eval_cols <- setdiff(colnames(eval), non_feature_cols)

  # identify correct gene columns
  common_genes <- intersect(training_cols, eval_cols)
  
  # old version, using padding
  #missing_genes <- setdiff(training_cols, eval_cols)
  #for (g in missing_genes) {
  #  eval[[g]] <- NA
  #}
  
  # subset correct columns
  eval <- eval[, c("patientID", "time", "event", "treatment", common_genes), drop=FALSE]
  
  # create RST30 column, in case it's needed
  tau <- 360
  eval$rst_30 <- mapply(
    function(time, event) {
      if (is.na(time) || is.na(event)) {
        return(NA)
      } else if (event == 1 && time <= tau) {
        return(time)
      } else {
        return(min(time, tau))
      }
    },
    eval$time, eval$event
  )
  tau <- 5
  eval$rst_5 <- mapply(
    function(time, event) {
      if (is.na(time) || is.na(event)) {
        return(NA)
      } else if (event == 1 && time <= tau) {
        return(time)
      } else {
        return(min(time, tau))
      }
    },
    eval$time, eval$event
  )
  
  return (eval)
}
