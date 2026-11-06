# Metrics from the raw estimates saved by run.R (no refitting).
# Usage: Rscript analyse.R
# Writes to RES_DIR:
#   metrics.rds   one row per replicate x version x cluster
#   summary.csv   means by scenario x cluster x version, with Monte Carlo errors

if (!exists("RES_DIR")) RES_DIR <- "results"
hosp <- readRDS(file.path(RES_DIR, "hosp_all.rds"))
clus <- readRDS(file.path(RES_DIR, "clus_all.rds"))
z95 <- qnorm(0.975)
KEYS <- c("J", "sigma_b", "sep", "rep", "version", "cluster")

# ---- hospitals: intervals as in the paper (REsim mean +- 1.96 sd) ----
h <- hosp
lo <- h$b_resim_mean - z95 * h$b_resim_sd
hi <- h$b_resim_mean + z95 * h$b_resim_sd
flag <- lo > 0 | hi < 0                          # hospital flagged as in Fig. 2
h$cov_b      <- lo <= h$b_true & h$b_true <= hi  # does the interval contain the true b?
h$width_b    <- hi - lo
h$flag_rate  <- flag
h$wrong_flag <- flag & sign(h$b_resim_mean) != sign(h$b_true)
h$missed     <- !flag & abs(h$b_true) >= h$sigma_b   # true effect >= 1 sigma_b, not flagged
hcols <- c("cov_b", "width_b", "flag_rate", "wrong_flag", "missed")
hm <- aggregate(h[, hcols], h[, KEYS], mean)     # shares per replicate

# ---- clusters: fixed effects and diagnostics ----
cl <- clus
cl$beta_cov <- abs(cl$beta_hat - cl$beta_true) <= z95 * cl$beta_se
cl$beta_err <- cl$beta_hat - cl$beta_true
cl$sb_zero  <- cl$sigma_b_hat < 1e-4
metrics <- merge(cl, hm, by = KEYS)
saveRDS(metrics, file.path(RES_DIR, "metrics.rds"))

# ---- summary by scenario x cluster x version ----
SCEN <- c("J", "sigma_b", "sep", "cluster", "version")
cov_cols <- c("beta_cov", "cov_b")
mean_cols <- c(cov_cols, "beta_err", "beta_se", "width_b", "flag_rate", "wrong_flag",
               "missed", "sb_zero", "contamination", "ari")
summ <- aggregate(metrics[, mean_cols], metrics[, SCEN], mean)
# Monte Carlo error of a coverage: sd across replicates / sqrt(number of replicates).
# The replicate is the unit: the hospitals of one replicate share the same
# estimated sigma_b and are correlated. For beta this equals sqrt(p(1-p)/R).
se <- aggregate(metrics[, cov_cols], metrics[, SCEN], function(x) sd(x) / sqrt(length(x)))
names(se)[-(1:5)] <- paste0(cov_cols, "_se")
n_reps <- aggregate(cbind(n_reps = metrics$beta_cov), metrics[, SCEN], length)
summ <- merge(merge(summ, se), n_reps)
summ <- summ[order(-summ$sep, summ$sigma_b, summ$J, summ$cluster, summ$version), ]
write.csv(summ, file.path(RES_DIR, "summary.csv"), row.names = FALSE)

print_cov <- function(col, title) {
  cat("\n==", title, "==\n")
  s <- summ
  s$txt <- sprintf("%.2f (%.3f)", s[[col]], s[[paste0(col, "_se")]])
  print(reshape(s[, c("J", "sigma_b", "sep", "cluster", "version", "txt")],
                idvar = c("J", "sigma_b", "sep", "cluster"), timevar = "version",
                direction = "wide"), row.names = FALSE)
}
print_cov("beta_cov", "Coverage of beta (Monte Carlo error)")
print_cov("cov_b", "Coverage of b, REsim intervals (Monte Carlo error)")
