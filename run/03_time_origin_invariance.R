# ---------------------------------------------------------------------------
# 05_nuisance_adjusted_prototype.R
#
# Small prototype of the nuisance-orthogonalized slope diagnostic, on two
# conditions only: intercept-only and slope-only.
#
# The point is NOT to adopt a new method now. It is to find out whether the
# centred-time diagnostic (Option A) is a bare workaround or whether it gives
# substantively the same answer as a nuisance-adjusted score (Option B).
#
# PASS CRITERION (fixed in advance)
#   The adjusted slope instability statistic computed on raw t = (0,1,2) and on
#   centred t = (-1,0,1) must agree to numerical tolerance. This is an algebraic
#   property of the projection, so a failure indicts the implementation first.
#
# ALSO EXPECTED
#   Option A (centred time, naive parm) should agree CLOSELY but NOT exactly
#   with Option B, because under centred time the two score blocks are nearly
#   but not exactly orthogonal (empirical correlation ~0.02-0.03, not 0).
#
# Run from the project root:
#     Rscript run/03_time_origin_invariance.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("R/era_gee_fixed.R")
source("R/nuisance_adjusted_slope.R")
source("R/synthetic_data.R")

line <- function(ch = "-") cat(strrep(ch, 78), "\n")

B1     <- REF_TRUTH$b1
B2_OTH <- REF_TRUTH$b2_others
B2_REC <- REF_TRUTH$b2_recent
B1_REC <- B1 - c(0.40, 0.30, 0.25, 0.35)

CONDS <- list(`slope-only`     = list(b1_recent = B1,     b2_recent = B2_REC),
              `intercept-only` = list(b1_recent = B1_REC, b2_recent = B2_OTH))

MINSIZE_ROW <- 450L; MAXDEPTH <- 3L

centre_time <- function(d) { d$X_time <- d$X_time - mean(sort(unique(d$X_time))); d }

fit_tree <- function(d, YN, ZN, Q, fitf, parm = NULL) {
  fml <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                           paste(ZN, collapse = " + ")))
  mob(fml, data = d, fit = fitf, cluster = d$id,
      control = mob_control(verbose = FALSE, maxdepth = MAXDEPTH,
                            minsize = MINSIZE_ROW, parm = parm))
}
root_test <- function(tr) {
  ti <- nodeapply(tr, 1, function(n) info_node(n))[[1]]
  sv <- if (width(tr) == 1L) "(none)"
        else names(tr$data)[nodeapply(tr, 1, function(n) split_node(n)$varid)[[1]]]
  list(split = sv, nnodes = width(tr),
       stat = ti$test["statistic", "IMM"], p = ti$test["p.value", "IMM"],
       all_stat = ti$test["statistic", ], all_p = ti$test["p.value", ])
}

for (cn in names(CONDS)) {
  cc  <- CONDS[[cn]]
  ref <- make_reference_data(b1_recent = cc$b1_recent, b2_recent = cc$b2_recent)
  dat <- ref$dat; YN <- ref$YN; ZN <- ref$ZN; Q <- ref$Q
  fam <- rep("gaussian", Q)
  datc <- centre_time(dat)

  line("="); cat(sprintf("CONDITION: %s\n", toupper(cn))); line("=")

  # ---------------------------------------------------------------------
  # STEP 1. Algebraic invariance of the adjusted score itself (root node)
  # ---------------------------------------------------------------------
  cat("\n1. adjusted score, raw vs centred time (root node)\n")
  m_raw <- ERA_GEE_fixed(dat[, YN],  dat$X_time,  fam, "independence", y_scale = "raw")
  m_cen <- ERA_GEE_fixed(datc[, YN], datc$X_time, fam, "independence", y_scale = "raw")

  cat("   raw    B:\n"); print(round(m_raw$B, 4))
  cat("   centred B (intercept = level at the mean occasion; slope unchanged):\n")
  print(round(m_cen$B, 4))
  cat(sprintf("   max |slope_raw - slope_centred| = %.3e  (slopes are origin-free)\n",
              max(abs(m_raw$B["slope", ] - m_cen$B["slope", ]))))

  a_raw <- nuisance_adjust_slope(m_raw$estfun, dat$id,  Q)
  a_cen <- nuisance_adjust_slope(m_cen$estfun, datc$id, Q)
  cat(sprintf("   max |U_S.I(raw) - U_S.I(centred)|, wave level   = %.3e\n",
              max(abs(a_raw - a_cen))))
  U_raw <- rowsum(a_raw, dat$id,  reorder = FALSE)
  U_cen <- rowsum(a_cen, datc$id, reorder = FALSE)
  cat(sprintf("   max |U_S.I(raw) - U_S.I(centred)|, subject level = %.3e\n",
              max(abs(U_raw - U_cen))))
  cat(sprintf("   raw:     cor(U_I, U_S) matched-outcome = %s\n",
      paste(sprintf("%.3f", diag(cor(rowsum(m_raw$estfun, dat$id, reorder = FALSE))[
              1:Q, (Q+1):(2*Q)])), collapse = " ")))
  cat(sprintf("   centred: cor(U_I, U_S) matched-outcome = %s\n",
      paste(sprintf("%.3f", diag(cor(rowsum(m_cen$estfun, datc$id, reorder = FALSE))[
              1:Q, (Q+1):(2*Q)])), collapse = " ")))

  # ---------------------------------------------------------------------
  # STEP 2. Trees
  # ---------------------------------------------------------------------
  cat("\n2. slope-focused trees\n")
  f_adj_raw <- era_gee_mob_fit_slopeadj(dat,  YN, "X_time", fam)
  f_adj_cen <- era_gee_mob_fit_slopeadj(datc, YN, "X_time", fam)
  f_nai_raw <- era_gee_mob_fit(dat,  YN, "X_time", fam, "independence", y_scale = "raw")
  f_nai_cen <- era_gee_mob_fit(datc, YN, "X_time", fam, "independence", y_scale = "raw")

  runs <- list(
    `B: adjusted,  raw t`     = list(d = dat,  f = f_adj_raw, parm = NULL),
    `B: adjusted,  centred t` = list(d = datc, f = f_adj_cen, parm = NULL),
    `A: naive parm, centred t`= list(d = datc, f = f_nai_cen, parm = (Q+1):(2*Q)),
    `   naive parm, raw t`    = list(d = dat,  f = f_nai_raw, parm = (Q+1):(2*Q))
  )
  out <- list()
  for (rn in names(runs)) {
    r  <- runs[[rn]]
    tr <- fit_tree(r$d, YN, ZN, Q, r$f, r$parm)
    rt <- root_test(tr)
    out[[rn]] <- rt
    cat(sprintf("   %-26s -> %-8s nodes=%d  IMM stat=%9.3f  p=%.4g\n",
                rn, rt$split, rt$nnodes, rt$stat, rt$p))
  }

  cat("\n3. PASS CRITERION\n")
  d_inv <- abs(out[["B: adjusted,  raw t"]]$stat - out[["B: adjusted,  centred t"]]$stat)
  cat(sprintf("   |stat_adj(raw) - stat_adj(centred)| = %.3e   %s\n", d_inv,
              if (d_inv < 1e-6) "PASS -- time-origin invariant" else "FAIL -- suspect the implementation"))
  cat(sprintf("   all-Z max |stat difference|          = %.3e\n",
              max(abs(out[["B: adjusted,  raw t"]]$all_stat -
                      out[["B: adjusted,  centred t"]]$all_stat))))

  cat("\n4. OPTION A vs OPTION B\n")
  dAB <- abs(out[["A: naive parm, centred t"]]$stat - out[["B: adjusted,  raw t"]]$stat)
  cat(sprintf("   A (centred, naive) IMM stat = %9.3f\n", out[["A: naive parm, centred t"]]$stat))
  cat(sprintf("   B (adjusted)       IMM stat = %9.3f\n", out[["B: adjusted,  raw t"]]$stat))
  cat(sprintf("   difference = %.3f  (%.2f%% of B)\n", dAB, 100 * dAB / out[["B: adjusted,  raw t"]]$stat))
  cat(sprintf("   same substantive conclusion (split variable): %s\n",
              identical(out[["A: naive parm, centred t"]]$split, out[["B: adjusted,  raw t"]]$split)))
  cat(sprintf("   contaminated raw-time naive parm IMM stat = %9.3f  -> %s\n\n",
              out[["   naive parm, raw t"]]$stat, out[["   naive parm, raw t"]]$split))
}
