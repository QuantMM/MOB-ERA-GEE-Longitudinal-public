# ---------------------------------------------------------------------------
# run_all.R
#
# Runs the whole validation suite from the repository root.
#
#   Rscript run_all.R           everything except the 1,000-replication study
#   Rscript run_all.R --full    including it (about 84 minutes on 32 cores)
#
# No data is required: everything is generated from R/synthetic_data.R with a
# fixed seed.
# ---------------------------------------------------------------------------

FULL  <- if (exists("FULL")) FULL else ("--full" %in% commandArgs(TRUE))
CORES <- max(1L, min(32L, parallel::detectCores(logical = TRUE) - 2L))

need <- c("partykit", "MASS", "matrixcalc")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss))
  stop("missing packages: ", paste(miss, collapse = ", "),
       "\n  install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")

if (!dir.exists("results")) dir.create("results")

scripts <- c("run/00_regression_test.R",
             "run/01_worked_example.R",
             "run/02_cluster_equivalence.R",
             "run/03_time_origin_invariance.R",
             "run/04_diagnostics_table.R")

hr <- function(ch = "=") cat(strrep(ch, 78), "\n")

for (s in scripts) {
  hr(); cat("RUNNING  ", s, "\n"); hr()
  t0 <- Sys.time()
  ok <- tryCatch({ source(s, local = new.env(), echo = FALSE); TRUE },
                 error = function(e) { cat("ERROR: ", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("\n[%s] %s  (%.1f s)\n\n", if (ok) "ok" else "FAILED", s,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

if (FULL) {
  hr(); cat("RUNNING   run/05_repeated_simulation.R  (1000 reps, ", CORES, " cores)\n", sep = "")
  hr()
  system2(file.path(R.home("bin"), "Rscript"),
          c("run/05_repeated_simulation.R", "1000", CORES))
} else {
  hr()
  cat("SKIPPED   run/05_repeated_simulation.R  (use --full to include it)\n")
  cat("          summarising the shipped 1,000-replication result instead:\n")
  hr()
  source("R/synthetic_data.R")            # REF_TRUTH, used by the summariser
  o <- readRDS("results/simulation_R1000.rds")
  res <- o$res; N_REP <- o$n_rep
  source("run/05b_summarise_simulation.R", local = TRUE)
}

hr(); cat("done\n"); hr()
