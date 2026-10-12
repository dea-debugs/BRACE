# Stratified RMST evaluation of one revision run (paper Table: treatment stratified survival comparison).
# Usage: Rscript evaluate_run.R RUN_LABEL
# Reads revision_runs/RUN_LABEL.csv, writes revision_runs/RUN_LABEL_eval.json.
# RMST difference = Followed minus Not Followed at tau = 150 months (survRM2::rmst2), pooled and within each
# observed arm (no chemo: treatment == 0, chemo: treatment == 1).

library(survRM2)
library(jsonlite)

REVISION_DIR <- "/Users/sindypinero/Downloads/4_Deanna/BRACE/revision_runs/"
TAU <- 150

label <- commandArgs(trailingOnly = TRUE)[1]
d <- read.csv(paste0(REVISION_DIR, label, ".csv"))
d$followed <- as.numeric(d$treatment == d$recommendation)

rmst_diff <- function(x) {
  n_f <- sum(x$followed == 1); n_nf <- sum(x$followed == 0)
  out <- list(n = nrow(x), n_followed = n_f, n_not_followed = n_nf)
  if (n_f < 2 || n_nf < 2) return (c(out, list(diff = NA, lower = NA, upper = NA, p = NA, note = "group too small")))
  r <- tryCatch(rmst2(time = x$time, status = x$event, arm = x$followed, tau = TAU), error = function(e) e)
  if (inherits(r, "error")) return (c(out, list(diff = NA, lower = NA, upper = NA, p = NA, note = conditionMessage(r))))
  u <- r$unadjusted.result
  c(out, list(diff = u[1, 1], lower = u[1, 2], upper = u[1, 3], p = u[1, 4],
              rmst_followed = r$RMST.arm1$rmst[["Est."]], rmst_not_followed = r$RMST.arm0$rmst[["Est."]]))
}

res <- list(
  run_label = label,
  tau = TAU,
  n_rec_chemo = sum(d$recommendation == 1),
  n_rec_no_chemo = sum(d$recommendation == 0),
  crosstab_observed_by_recommended = list(
    t0_r0 = sum(d$treatment == 0 & d$recommendation == 0), t0_r1 = sum(d$treatment == 0 & d$recommendation == 1),
    t1_r0 = sum(d$treatment == 1 & d$recommendation == 0), t1_r1 = sum(d$treatment == 1 & d$recommendation == 1)),
  pooled = rmst_diff(d),
  no_chemo_arm = rmst_diff(d[d$treatment == 0, ]),
  chemo_arm = rmst_diff(d[d$treatment == 1, ])
)

write_json(res, paste0(REVISION_DIR, label, "_eval.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
for (g in c("no_chemo_arm", "chemo_arm", "pooled")) {
  x <- res[[g]]
  cat(sprintf("%-13s n=%d (F %d / NF %d)  RMST diff %.2f (%.2f, %.2f)  p=%.4g\n",
              g, x$n, x$n_followed, x$n_not_followed, x$diff, x$lower, x$upper, x$p))
}
