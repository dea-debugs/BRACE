"""CausalPFN bridge called from R Code/hte.R via system("python run_causalpfn.py").

Reads (from the working directory):
  causalpfn_train.csv  columns Y, W, then feature columns
  causalpfn_test.csv   feature columns (same order as in train)
Writes:
  causalpfn_predictions.csv  one column, risk_score = estimated CATE per test row

Uses causalpfn 0.1.4 CATEEstimator with default hyperparameters. The package has
no seed argument, so seed 123 is set globally for random, numpy and torch.

torch, faiss and scikit-learn each ship their own libomp. On macOS arm64 running
them multithreaded in one process segfaults, so OpenMP is limited to one thread
(set before those libraries load). The process ends with os._exit(0) after the
predictions are written so a hang in library teardown can never block R's
system() call.
"""
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"

import random
import sys

import numpy as np
import pandas as pd
import torch
from causalpfn import CATEEstimator

SEED = 123


def main():
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)

    train = pd.read_csv("causalpfn_train.csv")
    test = pd.read_csv("causalpfn_test.csv")

    y = train["Y"].to_numpy(dtype=np.float32)
    t = train["W"].to_numpy(dtype=np.float32)
    X_train_df = train.drop(columns=["Y", "W"])

    if list(X_train_df.columns) != list(test.columns):
        sys.exit("run_causalpfn.py: feature columns in train and test do not match")

    X_train = X_train_df.to_numpy(dtype=np.float32)
    X_test = test.to_numpy(dtype=np.float32)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    estimator = CATEEstimator(device=device)
    estimator.fit(X_train, t, y)
    cate = estimator.estimate_cate(X_test)

    pd.DataFrame({"risk_score": cate}).to_csv("causalpfn_predictions.csv", index=False)
    sys.stdout.flush()
    sys.stderr.flush()


if __name__ == "__main__":
    main()
    os._exit(0)
