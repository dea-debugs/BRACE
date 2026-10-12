# REVISION EDIT: full paths to the brace env Python and revision_runs/run_causalpfn.py
CAUSALPFN_PYTHON <- "/Users/sindypinero/miniconda3/envs/brace/bin/python"
CAUSALPFN_SCRIPT <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/run_causalpfn.py"
CAUSALPFN_IO_DIR <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/causalpfn_io"

hte_causal_pfn <- function(df_train, X_train, X_test, outcome = "rst_30") {
  # REVISION EDIT: exchange files live in CAUSALPFN_IO_DIR (run_causalpfn.py reads/writes its working directory)
  # one subfolder per run label so runs can execute in parallel
  io_dir <- if (exists("RUN_LABEL")) file.path(CAUSALPFN_IO_DIR, RUN_LABEL) else CAUSALPFN_IO_DIR
  dir.create(io_dir, showWarnings = FALSE, recursive = TRUE)
  train_file <- file.path(io_dir, "causalpfn_train.csv")
  test_file  <- file.path(io_dir, "causalpfn_test.csv")
  pred_file  <- file.path(io_dir, "causalpfn_predictions.csv")
unlink(pred_file)

  # Save datasets
  write.csv(
    cbind(
      # EDIT Y = df_train$rst_5
      Y = df_train[[outcome]]
      , W = df_train$treatment
      , X_train
    )
    , train_file
    , row.names = FALSE
  )

  write.csv(
    cbind(
      X_test
    )
    , test_file
    , row.names = FALSE
  )

  # Call CausalPFN Python script (no native R implementation)
  cmd <- paste("cd", shQuote(io_dir), "&&", shQuote(CAUSALPFN_PYTHON), shQuote(CAUSALPFN_SCRIPT))
t0 <- Sys.time()
  status <- system(cmd)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("CausalPFN call: %.1f s (exit status %d)", elapsed, status))
  if (status != 0 || !file.exists(pred_file)) stop("run_causalpfn.py failed")

  preds <- read.csv(pred_file)
  est_effect <- preds$risk_score
  attr(est_effect, "causalpfn_seconds") <- elapsed

  return (est_effect)
}

hte_causal_tree <- function(y_train, X_train, X_test) {
  min_node <- 30
  
  # Define Causal Tree
  df_ct <- data.frame(
    outcome = y_train$rst_30,     # numeric outcome
    treatment = y_train$treatment,
    X_train
  )
  
  ct <- causalTree(
    formula = outcome ~ .,
    data = df_ct,
    treatment = df_ct$treatment,
    split.Rule = "CT",
    split.Honest = TRUE,
    cv.option = "CT",
    split.Bucket = FALSE,
    xval = 5,
    minsize = min_node
  )
  
  # Remove from memory
  rm(df_ct)
  gc()
  
  # Predict survival using test set (estimated treatment effect, not survival time - i.e., interpret as risk score)
  est_effect <- as.numeric(predict(ct, newdata = X_test))
  
  return (est_effect)
}

hte_causal_forest <- function(y_train, X_train, X_test) {
  min_node <- 30
  num_trees <- 6000
  
  # Define Causal Forest
  W <- y_train$treatment
  Y.time <- y_train$time
  Y.event <- y_train$event
  
  sf <- causal_survival_forest(
    X = X_train,
    Y = Y.time,
    W = W,
    D = Y.event,
    # EDIT horizon = 5*365,
    horizon = 360,
    num.trees = num_trees,
    sample.fraction = 0.4,
    min.node.size = min_node
  )
  
  # Predict survival using test set
  est_effect <- as.numeric(predict(sf, newdata = X_test)$predictions)
  
  return (est_effect)
}