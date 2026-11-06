# One replicate: the four versions of the 2 x 2 design on the same data.
#
#                     sigma_b known   sigma_b estimated
#   true clusters          A                B
#   CEM clusters           C                D   (method of the paper)
#
# Returns raw estimates only; metrics are computed in analyse.R.
#   hosp: one row per version x cluster x hospital (REsim mean and sd, true b)
#   clus: one row per version x cluster (fixed effects, sigma_b, assignment)

suppressPackageStartupMessages(library(mclust))   # adjustedRandIndex
source("sim.R")
source("cem.R")

# raw estimates of fit f, compared with true cluster c
extract <- function(f, c, z_used, tr) {
  bt <- tr$b_true[, c]
  lev <- as.integer(names(f$b))
  rs <- resim_b(f$model, names(f$b))
  hosp <- data.frame(cluster = c, hosp = lev,
                     b_resim_mean = rs$mean, b_resim_sd = rs$sd,
                     b_true = bt[lev], row.names = NULL)
  clus <- data.frame(cluster = c,
                     alpha_hat = f$beta[1], alpha_se = f$beta_se[1],
                     beta_hat = f$beta[2], beta_se = f$beta_se[2],
                     beta_true = tr$truth$beta[c],
                     sigma_b_hat = f$sigma_b,
                     n_hosp = length(lev), b_missing = length(bt) - length(lev),
                     n_assigned = sum(z_used == c),
                     contamination = mean(tr$data$z_true[z_used == c] != c),
                     row.names = NULL)
  list(hosp = hosp, clus = clus)
}

run_one <- function(J = 30, sigma_b = 0.5, sep = 3.5, seed = 1) {
  tr <- sim_data(J = J, sigma_b = sigma_b, sep = sep, seed = seed)
  d <- tr$data; z_true <- d$z_true
  hosp <- list(); clus <- list()
  add <- function(v, fits, z, ari) {
    for (c in 1:2) {
      e <- extract(fits[[c]], c, z, tr)
      e$hosp$version <- v; e$clus$version <- v; e$clus$ari <- ari
      hosp[[length(hosp) + 1]] <<- e$hosp; clus[[length(clus) + 1]] <<- e$clus
    }
  }

  # A and B: true clusters
  for (v in c("A", "B")) {
    s <- if (v == "A") sigma_b else NULL
    add(v, lapply(1:2, function(c) fit_cluster(d[z_true == c, ], s)), z_true, 1)
  }
  # C and D: clusters estimated by the CEM
  for (v in c("C", "D")) {
    s <- if (v == "C") sigma_b else NULL
    fit <- tryCatch(cem(d, 2, s), error = function(e) NULL)
    if (is.null(fit)) next
    add(v, fit$fits, fit$z, adjustedRandIndex(fit$z, z_true))
  }
  info <- data.frame(J = J, sigma_b = sigma_b, sep = sep, seed = seed)
  list(hosp = cbind(do.call(rbind, hosp), info), clus = cbind(do.call(rbind, clus), info))
}
