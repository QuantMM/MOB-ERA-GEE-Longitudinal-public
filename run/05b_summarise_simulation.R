# ---------------------------------------------------------------------------
# 07b_summarise_simulation.R
#
# Summarises 05_repeated_simulation.R into four clearly separated blocks:
#
#   A. ROOT RECOVERY          correct variable / correct root category partition
#   B. TREE COMPLEXITY        over-splitting below a correct root, node counts,
#                             and the depth-2 split variables (sanity check only)
#   C. SLOPE DIAGNOSTIC       specificity: pure null vs the intercept-only
#      SPECIFICITY            nuisance condition
#   D. PARAMETER RECOVERY     bias and ABSOLUTE RMSE, under two separate
#                             conditionings -- correct root cut, and exact
#                             two-node tree
#
# Coefficients are refitted on the ROOT's two branches in 05_repeated_simulation.R,
# so block D is available whenever the root cut is correct, including
# replications that over-split below it.
#
# Sourced at the end of that script; can also be run standalone:
#     o <- readRDS("results/simulation_R1000.rds"); res <- o$res; N_REP <- o$n_rep
#     source("run/05b_summarise_simulation.R")
# ---------------------------------------------------------------------------

if (!exists("res")) stop("no `res` in scope")
if (!exists("N_REP")) N_REP <- max(res$rep)

line <- function(ch = "-") cat(strrep(ch, 78), "\n")
rate <- function(x) 100 * mean(x, na.rm = TRUE)
mcse <- function(x) { p <- mean(x, na.rm = TRUE); 100 * sqrt(p * (1 - p) / sum(!is.na(x))) }
fmt  <- function(x) sprintf("%5.1f%% (MCSE %.1f)", rate(x), mcse(x))

ZN_ALL <- c("IMM", "SEX", "ETHN", "EDU", "INCNEED", "HOMEOWN", "URBAN")
Q <- 4L
YLAB <- c("walk", "lsport", "msport", "ssport")

POWER_CONDS <- c("slope-only", "intercept+slope")   # slopes truly differ
NUIS_COND   <- "intercept-only"                     # slopes homogeneous; intercepts differ

# ===========================================================================
# A. ROOT RECOVERY
# ===========================================================================
line("="); cat("A. ROOT RECOVERY\n"); line("=")
cat(sprintf("\n%d replications. Percentages carry Monte Carlo standard errors, not CIs.\n", N_REP))

cat("\nA1. Global null -- false split rate (nominal alpha = 0.05, Bonferroni over 7 variables)\n\n")
nul <- subset(res, condition == "null" & !error)
for (dn in c("joint", "slope", "intercept")) {
  d <- subset(nul, diagnostic == dn); if (!nrow(d)) next
  any_split <- d$split != "(none)"
  cat(sprintf("  %-10s P(any split) = %s\n", dn, fmt(any_split)))
  tb <- table(factor(d$split[any_split], levels = ZN_ALL))
  cat(sprintf("             by variable: %s\n",
              paste(sprintf("%s=%.1f%%", names(tb), 100 * tb / nrow(d)), collapse = "  ")))
}

cat("\nA2. Root variable and root category partition, non-null conditions\n")
cat("    The intercept-only / SLOPE cell is deliberately absent: the slopes are\n")
cat("    homogeneous there by construction, so it is a specificity result, not a\n")
cat("    recovery result. It is reported in block C.\n")
cat("    NOTE: power here is CEILING power at a deliberately large planted signal;\n")
cat("    it is not a general claim about the method's sensitivity.\n\n")
pw <- subset(res, condition %in% POWER_CONDS & !error)
aggset <- subset(res, !error & ((condition %in% POWER_CONDS) |
                 (condition == NUIS_COND & diagnostic == "joint")))
aggA <- do.call(rbind, lapply(
  split(aggset, list(aggset$condition, aggset$diagnostic), drop = TRUE), function(d)
  data.frame(condition = d$condition[1], diagnostic = d$diagnostic[1], n = nrow(d),
             root_var_IMM   = round(rate(d$root_var_correct), 1),
             wrong_variable = round(rate(d$split != "IMM" & d$split != "(none)"), 1),
             no_split       = round(rate(d$split == "(none)"), 1),
             root_cut_exact = round(rate(d$root_cut_correct), 1))))
print(aggA[order(aggA$condition, aggA$diagnostic), ], row.names = FALSE)
cat("\n  root_cut_exact = root partition is exactly {recent} | {established, non}\n")

# ===========================================================================
# B. TREE COMPLEXITY
# ===========================================================================
line("="); cat("B. TREE COMPLEXITY\n"); line("=")
cat("\nB1. Terminal-node counts\n\n")
for (cnd in c(POWER_CONDS, NUIS_COND, "null")) {
  d <- subset(res, condition == cnd & !error)
  for (dn in unique(d$diagnostic)) {
    dd <- subset(d, diagnostic == dn)
    tb <- table(dd$n_nodes)
    cat(sprintf("  %-16s %-10s %s\n", cnd, dn,
                paste(sprintf("%s:%.1f%%", names(tb), 100 * tb / nrow(dd)), collapse = "  ")))
  }
}

cat("\nB2. Over-splitting CONDITIONAL on a correct root cut\n\n")
aggB <- do.call(rbind, lapply(split(pw, list(pw$condition, pw$diagnostic), drop = TRUE), function(d) {
  s <- subset(d, root_cut_correct)
  data.frame(condition = d$condition[1], diagnostic = d$diagnostic[1],
             n_correct_root = nrow(s),
             exact_two_node = round(rate(s$exact_two_node_tree), 1),
             oversplit      = round(rate(s$oversplit_after_correct_root), 1),
             mean_nodes     = round(mean(s$n_nodes), 3))
}))
print(aggB[order(aggB$condition, aggB$diagnostic), ], row.names = FALSE)

cat("\nB3. Depth-2 split variables, among replications that over-split a correct root\n")
cat("    (supplementary sanity check, not a primary result)\n\n")
os <- subset(pw, oversplit_after_correct_root)
if (nrow(os)) {
  d2 <- c(os$depth2_a, os$depth2_b)
  d2 <- d2[d2 != "(none)" & !is.na(d2)]
  tb <- sort(table(factor(d2, levels = ZN_ALL)), decreasing = TRUE)
  cat(sprintf("  %d secondary splits over %d over-splitting replications\n", length(d2), nrow(os)))
  print(tb)
  top <- 100 * max(tb) / sum(tb)
  nz  <- sum(tb > 0)
  cat(sprintf("\n  %d of 7 variables appear; the most frequent accounts for %.1f%%\n", nz, top))
  cat(sprintf("  expected share if secondary splits were uniform over 7 variables: %.1f%%\n",
              100 / 7))
  if (nz >= 4 && top < 40)
    cat("  => scattered across several null variables, consistent with stochastic\n     over-splitting rather than a reproducible secondary signal.\n")
  else
    cat("  => concentrated. Check whether a specific variable is being picked up\n     repeatedly before calling this stochastic over-splitting.\n")
} else cat("  none\n")

# ===========================================================================
# C. SLOPE DIAGNOSTIC SPECIFICITY
# ===========================================================================
line("="); cat("C. SLOPE DIAGNOSTIC SPECIFICITY\n"); line("=")
cat("\nUnder intercept-only truth the slopes are homogeneous by construction, so\n")
cat("for the SLOPE diagnostic that condition is itself a null: any split is a\n")
cat("false positive. This is the check that intercept heterogeneity does not\n")
cat("contaminate the nuisance-adjusted slope test.\n\n")
sn <- subset(res, condition == NUIS_COND & diagnostic == "slope" & !error)
s0 <- subset(res, condition == "null"    & diagnostic == "slope" & !error)
cat(sprintf("  pure null                      P(any false split) = %s\n",
            fmt(s0$split != "(none)")))
cat(sprintf("  intercept-only (nuisance)      P(any false split) = %s\n",
            fmt(sn$split != "(none)")))
cat(sprintf("     of which on IMM: %.1f%%   on some other variable: %.1f%%\n",
            rate(sn$split == "IMM"), rate(sn$split != "IMM" & sn$split != "(none)")))
cat("\n  For contrast, the joint diagnostic splits in this condition by design,\n")
cat("  because the trajectories genuinely differ in baseline:\n")
jn <- subset(res, condition == NUIS_COND & diagnostic == "joint" & !error)
cat(sprintf("  joint, intercept-only          P(split on IMM)    = %s\n",
            fmt(jn$split == "IMM")))

# ===========================================================================
# D. PARAMETER RECOVERY
# ===========================================================================
line("="); cat("D. PARAMETER RECOVERY\n"); line("=")
cat("\nCoefficients are refitted on the ROOT's two branches, so they are defined\n")
cat("whenever the root cut is correct, over-splitting or not. Two conditionings\n")
cat("are reported separately.\n")
cat("\nABSOLUTE RMSE is the primary quantity. Relative RMSE is not reported as a\n")
cat("headline: the true slopes differ by two orders of magnitude across conditions\n")
cat("(-0.45 down to -0.02), so a ratio explodes wherever the truth is near zero.\n")
cat("\nThe four conditions share replication seeds (common random numbers) and differ\n")
cat("only by a deterministic shift of the 'recent' group's truth, so bias and RMSE\n")
cat("come out identical across conditions. Expected, not a bug.\n")

TRUTH <- list(
  `slope-only`      = list(b1 = REF_TRUTH$b1,
                           b2 = REF_TRUTH$b2_recent),
  `intercept-only`  = list(b1 = REF_TRUTH$b1 - c(0.40, 0.30, 0.25, 0.35),
                           b2 = REF_TRUTH$b2_others),
  `intercept+slope` = list(b1 = REF_TRUTH$b1 - c(0.40, 0.30, 0.25, 0.35),
                           b2 = REF_TRUTH$b2_recent)
)

recov <- function(d, cn, label) {
  if (!nrow(d)) return(invisible(NULL))
  cat(sprintf("\n-- %s | %s  (n = %d)\n", cn, label, nrow(d)))
  tb <- do.call(rbind, lapply(1:Q, function(q) {
    e <- function(col, truth) {
      v <- d[[col]]
      c(bias = mean(v - truth, na.rm = TRUE), rmse = sqrt(mean((v - truth)^2, na.rm = TRUE)))
    }
    r1 <- e(paste0("b1_rec_", q), TRUTH[[cn]]$b1[q])
    r2 <- e(paste0("b2_rec_", q), TRUTH[[cn]]$b2[q])
    r3 <- e(paste0("b1_oth_", q), REF_TRUTH$b1[q])
    r4 <- e(paste0("b2_oth_", q), REF_TRUTH$b2_others[q])
    data.frame(outcome = YLAB[q],
               true_b2_recent = TRUTH[[cn]]$b2[q],
               b1rec_bias = round(r1[1], 4), b1rec_rmse = round(r1[2], 4),
               b2rec_bias = round(r2[1], 4), b2rec_rmse = round(r2[2], 4),
               b1oth_bias = round(r3[1], 4), b1oth_rmse = round(r3[2], 4),
               b2oth_bias = round(r4[1], 4), b2oth_rmse = round(r4[2], 4))
  }))
  print(tb, row.names = FALSE)
}

for (cn in names(TRUTH)) {
  for (dn in c("joint", "slope")) {
    base <- subset(res, condition == cn & diagnostic == dn & !error)
    recov(subset(base, root_cut_correct),    cn, sprintf("%s | conditional on CORRECT ROOT CUT", dn))
    recov(subset(base, exact_two_node_tree), cn, sprintf("%s | conditional on EXACT TWO-NODE TREE", dn))
  }
}

line("="); cat("TIMING\n"); line("=")
cat(sprintf("mean %.1f s per tree | median %.1f | max %.1f\n",
            mean(res$secs), median(res$secs), max(res$secs)))
print(round(tapply(res$secs, list(res$condition, res$diagnostic), mean), 1))
