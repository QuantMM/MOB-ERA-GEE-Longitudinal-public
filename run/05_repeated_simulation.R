# ---------------------------------------------------------------------------
# 07_repeated_simulation.R
#
# Repeated simulation for the adopted MOB-ERA-GEE diagnostics.
#
# ADOPTED DEFINITIONS (all on raw t = 0,1,2, all with cluster = id):
#   joint      full 2Q model, instability test on all 2Q parameters
#   slope      full 2Q model, instability test on the Q nuisance-adjusted
#              slope scores  U_S.I = U_S - I_SI I_II^-1 U_I
#   intercept  full 2Q model, instability test on the Q intercepts (parm = 1:Q)
#
# The intercept diagnostic is a supplementary check, not a methodological
# contribution, so it is run only under the NULL condition, where its
# calibration matters. joint and slope are run in all four conditions.
#
# FOUR SEPARATED QUANTITIES
#   1. false split rate under the null
#        - P(any split)         -- should sit near alpha with Bonferroni
#        - P(split on each specific variable)
#   2. true partitioning-variable recovery (non-null)
#        - P(first split on IMM), plus the wrong-variable rate and the miss rate
#   3. correct category split recovery
#        - P(the binary partition is exactly {recent} | {non, established})
#   4. terminal-node coefficient bias / RMSE for beta_1q and beta_2q,
#      reported conditional on correct recovery
#   (plus the distribution of the number of terminal nodes, i.e. over-splitting)
#
# Usage:
#     Rscript run/05_repeated_simulation.R [n_rep] [n_cores]
# Defaults to a short pilot so the harness and the per-replication cost can be
# checked before the full run.
# ---------------------------------------------------------------------------

args    <- commandArgs(trailingOnly = TRUE)
N_REP   <- if (length(args) >= 1) as.integer(args[1]) else 20L
N_CORES <- if (length(args) >= 2) as.integer(args[2]) else 16L
SEED0   <- 700000L

N_PERSON    <- 3000L
MINSIZE_ROW <- 450L
MAXDEPTH    <- 3L

suppressPackageStartupMessages({library(partykit); library(MASS); library(parallel)})
source("R/synthetic_data.R")

B1     <- REF_TRUTH$b1
B2_OTH <- REF_TRUTH$b2_others
B2_REC <- REF_TRUTH$b2_recent
B1_REC <- B1 - c(0.40, 0.30, 0.25, 0.35)

CONDITIONS <- list(
  `slope-only`      = list(b1_recent = B1,     b2_recent = B2_REC, diags = c("joint","slope")),
  `intercept-only`  = list(b1_recent = B1_REC, b2_recent = B2_OTH, diags = c("joint","slope")),
  `intercept+slope` = list(b1_recent = B1_REC, b2_recent = B2_REC, diags = c("joint","slope")),
  `null`            = list(b1_recent = B1,     b2_recent = B2_OTH,
                           diags = c("joint","slope","intercept"))
)

# ---------------------------------------------------------------------------
# one replication x one condition -> one row per diagnostic
# ---------------------------------------------------------------------------
one_rep <- function(job) {
  cn <- job$cond; rep_i <- job$rep
  cc <- CONDITIONS[[cn]]
  ref <- make_reference_data(n_person = N_PERSON, seed = SEED0 + rep_i,
                             b1_recent = cc$b1_recent, b2_recent = cc$b2_recent)
  dat <- ref$dat; YN <- ref$YN; ZN <- ref$ZN; Q <- ref$Q
  fam <- rep("gaussian", Q)
  fml <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                           paste(ZN, collapse = " + ")))

  truth_b1_rec <- cc$b1_recent; truth_b2_rec <- cc$b2_recent

  out <- list()
  for (dn in cc$diags) {
    fitf <- if (dn == "slope")
      era_gee_mob_fit_slopeadj(dat, YN, "X_time", fam)
    else
      era_gee_mob_fit(dat, YN, "X_time", fam, "independence", y_scale = "raw")
    prm <- switch(dn, joint = 1:(2 * Q), intercept = 1:Q, slope = NULL)

    t0 <- Sys.time()
    tr <- try(mob(fml, data = dat, fit = fitf, cluster = dat$id,
                  control = mob_control(verbose = FALSE, maxdepth = MAXDEPTH,
                                        minsize = MINSIZE_ROW, parm = prm)),
              silent = TRUE)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (inherits(tr, "try-error")) {
      out[[dn]] <- data.frame(condition = cn, rep = rep_i, diagnostic = dn,
                              error = TRUE, split = NA_character_, n_nodes = NA_integer_,
                              root_var_correct = NA, root_cut_correct = NA,
                              exact_two_node_tree = NA, oversplit_after_correct_root = NA,
                              root_left = NA_character_, root_right = NA_character_,
                              depth2_a = NA_character_, depth2_b = NA_character_,
                              secs = el, stringsAsFactors = FALSE)
      next
    }

    nn   <- width(tr)
    root <- node_party(tr)
    sv   <- if (nn == 1L) "(none)" else names(tr$data)[split_node(root)$varid]

    # ------------------------------------------------------------------
    # ROOT PARTITION, stored as human-readable CATEGORY LABELS.
    # Numeric split indices are not archived on their own: if factor level
    # ordering ever changes they become uninterpretable. The label strings
    # plus the seed are enough to reconstruct root membership later, even
    # without the tree object.
    # ------------------------------------------------------------------
    root_left <- root_right <- NA_character_
    root_kid  <- NULL
    if (nn > 1L) {
      sp  <- split_node(root)
      zv  <- tr$data[[sp$varid]]
      if (!is.null(sp$index) && is.factor(zv)) {
        lv  <- levels(zv)
        idx <- sp$index
        root_left  <- paste(sort(lv[which(idx == 1L)]), collapse = "|")
        root_right <- paste(sort(lv[which(idx == 2L)]), collapse = "|")
        root_kid   <- idx[as.integer(dat[[sv]])]        # 1 / 2 per row
      } else if (!is.null(sp$breaks)) {                 # numeric splitter
        root_left  <- paste0("<=", signif(sp$breaks[1], 6))
        root_right <- paste0(">",  signif(sp$breaks[1], 6))
        root_kid   <- ifelse(dat[[sv]] <= sp$breaks[1], 1L, 2L)
      }
    }

    root_var_correct <- identical(sv, "IMM")
    root_cut_correct <- root_var_correct &&
      ((identical(root_left, "recent")  && identical(root_right, "established|non")) ||
       (identical(root_right, "recent") && identical(root_left,  "established|non")))
    exact_two_node_tree <- root_cut_correct && nn == 2L
    oversplit_after_correct_root <- root_cut_correct && nn > 2L

    # depth-2 split variables, one per root branch (supplementary sanity check)
    d2 <- c("(none)", "(none)")
    if (nn > 1L) {
      kd <- kids_node(root)
      d2 <- sapply(seq_along(kd), function(k) {
        s <- split_node(kd[[k]])
        if (is.null(s)) "(none)" else names(tr$data)[s$varid]
      })
      if (length(d2) < 2L) d2 <- c(d2, rep("(none)", 2L - length(d2)))
    }

    # ------------------------------------------------------------------
    # Coefficients are refitted on the ROOT's two branches, not on the
    # terminal nodes, so that root-level recovery can be assessed whenever
    # the root cut is correct -- including replications that over-split
    # below it.
    # ------------------------------------------------------------------
    b1e <- b2e <- rep(NA_real_, 2 * Q)
    if (isTRUE(root_cut_correct)) {
      rec_kid <- if (identical(root_left, "recent")) 1L else 2L
      for (k in 1:2) {
        d <- dat[root_kid == k, ]
        m <- ERA_GEE_fixed(d[, YN], d$X_time, fam, "independence", y_scale = "raw")
        if (k == rec_kid) { b1e[1:Q] <- m$B["intercept", ]; b2e[1:Q] <- m$B["slope", ] }
        else { b1e[(Q+1):(2*Q)] <- m$B["intercept", ]; b2e[(Q+1):(2*Q)] <- m$B["slope", ] }
      }
    }

    r <- data.frame(condition = cn, rep = rep_i, diagnostic = dn, error = FALSE,
                    split = sv, n_nodes = nn,
                    root_var_correct = root_var_correct,
                    root_cut_correct = root_cut_correct,
                    exact_two_node_tree = exact_two_node_tree,
                    oversplit_after_correct_root = oversplit_after_correct_root,
                    root_left = root_left, root_right = root_right,
                    depth2_a = d2[1], depth2_b = d2[2],
                    secs = el, stringsAsFactors = FALSE)
    est <- as.data.frame(t(c(b1e, b2e)))
    names(est) <- c(paste0("b1_rec_", 1:Q), paste0("b1_oth_", 1:Q),
                    paste0("b2_rec_", 1:Q), paste0("b2_oth_", 1:Q))
    out[[dn]] <- cbind(r, est)
  }
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
jobs <- do.call(c, lapply(names(CONDITIONS), function(cn)
  lapply(seq_len(N_REP), function(i) list(cond = cn, rep = i))))

cat(sprintf("repeated simulation: %d reps x %d conditions = %d jobs, %d cores\n",
            N_REP, length(CONDITIONS), length(jobs), N_CORES))
cat(sprintf("n_person = %d, maxdepth = %d, minsize = %d rows\n\n",
            N_PERSON, MAXDEPTH, MINSIZE_ROW))

t_start <- Sys.time()
cl <- makeCluster(min(N_CORES, length(jobs)))
on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(partykit); library(MASS)})
  source("R/era_gee_fixed.R")
  source("R/nuisance_adjusted_slope.R")
  source("R/synthetic_data.R")
  NULL
}))
clusterExport(cl, c("CONDITIONS", "N_PERSON", "MINSIZE_ROW", "MAXDEPTH", "SEED0",
                    "B1", "B2_OTH", "B2_REC", "B1_REC", "one_rep"))
res <- do.call(rbind, parLapplyLB(cl, jobs, one_rep))
stopCluster(cl)
elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

cat(sprintf("done in %.1f min | mean seconds per tree = %.1f | errors = %d\n\n",
            elapsed, mean(res$secs, na.rm = TRUE), sum(res$error)))

saveRDS(list(res = res, n_rep = N_REP, n_person = N_PERSON,
             maxdepth = MAXDEPTH, minsize = MINSIZE_ROW, elapsed_min = elapsed),
        sprintf("results/_simulation_R%d.rds", N_REP))

source("run/05b_summarise_simulation.R", local = TRUE)
