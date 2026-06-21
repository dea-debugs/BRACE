hte_causal_pfn <- function(df_train, X_train, X_test) {
  # Save datasets
  write.csv(
    cbind(
      # EDIT Y = df_train$rst_5
      Y = df_train$rst_30
      , W = df_train$treatment
      , X_train
    )
    , "causalpfn_train.csv"
    , row.names = FALSE
  )
  
  write.csv(
    cbind(
      X_test
    )
    , "causalpfn_test.csv"
    , row.names = FALSE
  )
  
  # Call CausalPFN Python script (no native R implementation)
  system("python run_causalpfn.py")
  preds <- read.csv("causalpfn_predictions.csv")
  est_effect <- preds$risk_score
  
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