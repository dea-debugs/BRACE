# grf causal survival forest for the semi synthetic benchmark (run F).
# Usage: Rscript csf_fit.R TRAIN_CSV TEST_CSV OUT_CSV SEED HORIZON NUM_THREADS
# TRAIN_CSV: U (observed time), D (event), W (treatment), then features; TEST_CSV: same features.
# Writes OUT_CSV with column tau_hat = estimated RMST(horizon) difference, treated minus control.

suppressMessages(library(grf))

a <- commandArgs(trailingOnly = TRUE)
train <- read.csv(a[1]); test <- read.csv(a[2])
seed <- as.integer(a[4]); horizon <- as.numeric(a[5]); threads <- as.integer(a[6])

X <- as.matrix(train[, setdiff(colnames(train), c("U", "D", "W"))])
Xt <- as.matrix(test[, colnames(X)])

csf <- causal_survival_forest(X, Y = train$U, W = train$W, D = train$D,
                              target = "RMST", horizon = horizon,
                              num.trees = 2000, seed = seed, num.threads = threads)
write.csv(data.frame(tau_hat = predict(csf, Xt)$predictions), a[3], row.names = FALSE)
