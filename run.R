# Full simulation: 8 scenarios x 100 replicates, versions A-D.
# Usage: Rscript run.R   (SMOKE=1 for a quick test in results_smoke/)
# Each replicate is saved in results/ as soon as it is done; replicates already
# there are not recomputed. analyse.R is run at the end.

source("study.R")
suppressPackageStartupMessages(library(parallel))

R_REPS    <- 100
N_WORKERS <- max(1, detectCores() - 1)
RES_DIR   <- "results"
grid <- expand.grid(J       = c(10, 30),     # hospitals
                    sigma_b = c(0.5, 1.5),   # similar / very different hospitals
                    sep     = c(3.5, 2),     # separated as in the application / overlapping
                    rep     = 1:R_REPS)
if (Sys.getenv("SMOKE") == "1") {
  grid <- expand.grid(J = c(10, 30), sigma_b = 0.5, sep = 3.5, rep = 1:2)
  RES_DIR <- "results_smoke"
}

dir.create(RES_DIR, showWarnings = FALSE)
run_row <- function(k) {
  g <- grid[k, ]
  f <- sprintf("%s/J%02d_sb%.1f_sep%.1f_rep%03d.rds", RES_DIR, g$J, g$sigma_b, g$sep, g$rep)
  if (file.exists(f)) return(readRDS(f))
  seed <- g$rep + 1000 * g$J + round(1e5 * g$sigma_b) + round(1e6 * g$sep)
  res <- tryCatch(run_one(g$J, g$sigma_b, g$sep, seed),
                  error = function(e) { message("error ", f, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) { res$hosp$rep <- g$rep; res$clus$rep <- g$rep; saveRDS(res, f) }
  res
}

t_start <- Sys.time()
out <- mclapply(seq_len(nrow(grid)), run_row, mc.cores = N_WORKERS)
out <- out[!sapply(out, is.null)]
saveRDS(do.call(rbind, lapply(out, `[[`, "hosp")), file.path(RES_DIR, "hosp_all.rds"))
saveRDS(do.call(rbind, lapply(out, `[[`, "clus")), file.path(RES_DIR, "clus_all.rds"))
cat(sprintf("Simulation: %.1f minutes\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

source("analyse.R")
