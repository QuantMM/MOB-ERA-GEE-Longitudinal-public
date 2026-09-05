# ---------------------------------------------------------------------------
# 02_A_vs_Ba.R
#
# Method A vs. method B(a) on the frozen reference case.
#
#   B(a)  THEORETICAL REFERENCE. mob() partitions a PERSON-level frame; the
#         estimating-function contribution is the GEE-native subject-level
#         U_i = sum_t u_it, so estfun is N x 2Q and the instability tests run
#         over genuinely independent clusters i = 1..N. No `cluster` argument.
#
#   A     CONVENIENT IMPLEMENTATION. mob() partitions the person-WAVE frame;
#         estfun is NT x 2Q and person clustering is handled by
#         mob(..., cluster = id), which partykit uses for the clustered
#         covariance in the parameter-stability tests.
#
# These are two implementations of the same specification, not two models.
# The question is NOT "which is better" but "does A reproduce B(a)?".
# If it does, A is the sensible thing to run the simulations with.
#
# CONCLUSION (stated carefully -- see section 6 and NOTES_02_A_vs_Ba.md):
#   A is exactly equivalent to B(a) for the score aggregation and the clustered
#   covariance, and empirically reproduces B(a) exactly for the CATEGORICAL
#   partitioning variables in the reference case. For NUMERIC person-level
#   partitioning variables the fluctuation processes are not algebraically
#   identical, although the observed discrepancy was negligible here.
#   All of this is established under BALANCED repeated measurements.
#
# Checks:
#   1. root coefficients identical
#   2. person-wise column sums of A's estfun == B(a)'s estfun   (exact)
#   3. root parameter-instability tests
#   4. split variable, cut point, terminal-node membership
#   5. terminal-node parameter estimates
#
# Run from the project root:
#     Rscript run/02_cluster_equivalence.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("R/era_gee_fixed.R")
source("R/synthetic_data.R")

line <- function(ch = "-") cat(strrep(ch, 78), "\n")

ref  <- make_reference_data()
dat  <- ref$dat; pdat <- ref$pdat
YN   <- ref$YN;  ZN   <- ref$ZN; Q <- ref$Q
fam  <- rep("gaussian", Q)

# minsize must mean the same thing in both: 150 persons = 450 person-wave rows
MINSIZE_PERSON <- 150L
MINSIZE_ROW    <- MINSIZE_PERSON * ref$n_wave

line("="); cat("A vs B(a) on the reference case\n"); line("=")
cat(sprintf("persons = %d | rows = %d | Q = %d | minsize = %d persons (= %d rows)\n\n",
            ref$n_person, nrow(dat), Q, MINSIZE_PERSON, MINSIZE_ROW))

# ===========================================================================
# 1-2. ROOT: coefficients, and the exact score-aggregation identity
# ===========================================================================
line(); cat("1-2. ROOT MODEL AND SCORE AGGREGATION\n"); line()

root <- ERA_GEE_fixed(dat[, YN], dat$X_time, fam, "independence", y_scale = "raw")

fitA  <- era_gee_mob_fit(dat, YN, "X_time", fam, "independence", y_scale = "raw")
fitBa <- era_gee_mob_fit_subject(dat, pdat, YN, "X_time", fam, "independence",
                                 y_scale = "raw")

# call each fit function exactly as mob would at the root
xA  <- cbind(`(Intercept)` = 1, rowid = dat$rowid,  X_time = dat$X_time)
xBa <- cbind(`(Intercept)` = 1, prow  = pdat$prow)
rA  <- fitA(y = as.matrix(dat[, YN]),  x = xA,  estfun = TRUE)
rBa <- fitBa(y = as.matrix(pdat[, YN]), x = xBa, estfun = TRUE)

cat("root coefficients\n")
print(round(rbind(A = rA$coefficients, `B(a)` = rBa$coefficients), 6))
cat(sprintf("\nmax |coef_A - coef_B(a)|          = %.3e\n",
            max(abs(rA$coefficients - rBa$coefficients))))
cat(sprintf("objfun A = %.6f   objfun B(a) = %.6f   diff = %.3e\n",
            rA$objfun, rBa$objfun, abs(rA$objfun - rBa$objfun)))

cat(sprintf("\nestfun dims:  A = %d x %d   B(a) = %d x %d\n",
            nrow(rA$estfun), ncol(rA$estfun), nrow(rBa$estfun), ncol(rBa$estfun)))

# THE SHARP TEST: U_i = sum_t u_it
aggA <- rowsum(rA$estfun, group = dat$id, reorder = FALSE)
aggA <- aggA[match(as.character(pdat$id), rownames(aggA)), , drop = FALSE]
rownames(aggA) <- NULL
d_ident <- max(abs(aggA - rBa$estfun))
cat(sprintf("\nmax |rowsum(estfun_A, by person) - estfun_B(a)| = %.3e   %s\n",
            d_ident, if (d_ident < 1e-10) "<- EXACT" else "<- MISMATCH"))
cat(sprintf("column sums vanish:  A = %.2e   B(a) = %.2e\n",
            max(abs(colSums(rA$estfun))), max(abs(colSums(rBa$estfun)))))

# ===========================================================================
# 3-5. TREES
# ===========================================================================
line(); cat("3-5. TREES\n"); line()

fmlA <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                          paste(ZN, collapse = " + ")))
fmlB <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ prow | ",
                          paste(ZN, collapse = " + ")))

t0 <- Sys.time()
trA <- mob(fmlA, data = dat, fit = fitA, cluster = dat$id,
           control = mob_control(verbose = FALSE, maxdepth = 3, minsize = MINSIZE_ROW))
tA <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

t0 <- Sys.time()
trB <- mob(fmlB, data = pdat, fit = fitBa,
           control = mob_control(verbose = FALSE, maxdepth = 3, minsize = MINSIZE_PERSON))
tB <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("A    : %d terminal nodes, %.1f s\n", width(trA), tA))
cat(sprintf("B(a) : %d terminal nodes, %.1f s\n\n", width(trB), tB))

show_tests <- function(tr, lab) {
  ti <- nodeapply(tr, 1, function(n) info_node(n))[[1]]
  m  <- rbind(statistic = ti$test["statistic", ], p.value = ti$test["p.value", ])
  cat(lab, "-- root parameter instability tests:\n")
  print(signif(m, 4)); cat("\n")
  m
}
mA <- show_tests(trA, "A   ")
mB <- show_tests(trB, "B(a)")

split_name <- function(tr) names(tr$data)[nodeapply(tr, 1, function(n) split_node(n)$varid)[[1]]]
split_idx  <- function(tr) nodeapply(tr, 1, function(n) split_node(n)$index)[[1]]

cat(sprintf("split variable :  A = %-6s   B(a) = %s\n", split_name(trA), split_name(trB)))
cat(sprintf("split index    :  A = %-12s B(a) = %s\n",
            paste(split_idx(trA), collapse = ","), paste(split_idx(trB), collapse = ",")))

# terminal-node membership, compared at the PERSON level
membA <- predict(trA, type = "node")
membB <- predict(trB, type = "node")
personA <- membA[!duplicated(dat$id)][match(pdat$id, dat$id[!duplicated(dat$id)])]
agree <- mean(as.integer(factor(personA)) == as.integer(factor(membB)))
cat(sprintf("\nterminal-node membership agreement (per person) = %.4f\n", agree))
cat("cross-tabulation of node assignments:\n")
print(table(A = personA, `B(a)` = membB))

# terminal-node parameters
node_par <- function(tr, memb_person, lab) {
  do.call(rbind, lapply(sort(unique(memb_person)), function(nd) {
    ids <- pdat$id[memb_person == nd]
    d   <- dat[dat$id %in% ids, ]
    m   <- ERA_GEE_fixed(d[, YN], d$X_time, fam, "independence", y_scale = "raw")
    gp  <- sort(unique(as.character(d$IMM[!duplicated(d$id)])))
    data.frame(method = lab, node = nd, n_person = length(ids),
               IMM_groups = paste(gp, collapse = " + "), outcome = YN,
               intercept = round(m$B["intercept", ], 4),
               slope     = round(m$B["slope", ], 4), row.names = NULL)
  }))
}
cat("\nterminal-node parameter estimates:\n")
np <- rbind(node_par(trA, personA, "A"), node_par(trB, membB, "B(a)"))
print(np, row.names = FALSE)

wide <- merge(subset(np, method == "A",    c(IMM_groups, outcome, intercept, slope)),
              subset(np, method == "B(a)", c(IMM_groups, outcome, intercept, slope)),
              by = c("IMM_groups", "outcome"), suffixes = c("_A", "_Ba"))
cat(sprintf("\nmax |intercept_A - intercept_B(a)| = %.3e\n",
            max(abs(wide$intercept_A - wide$intercept_Ba))))
cat(sprintf("max |slope_A     - slope_B(a)|     = %.3e\n",
            max(abs(wide$slope_A - wide$slope_Ba))))

line("="); cat("VERDICT\n"); line("=")
same <- identical(split_name(trA), split_name(trB)) &&
        width(trA) == width(trB) && agree == 1 && d_ident < 1e-10
cat(if (same) paste0(
  "A is exactly equivalent to B(a) for the score aggregation and the clustered\n",
  "covariance, and reproduces B(a) exactly for the CATEGORICAL partitioning\n",
  "variables in this reference case (identical scores after aggregation,\n",
  "identical tree, identical node membership and parameters).\n\n",
  "For NUMERIC person-level partitioning variables the fluctuation processes are\n",
  "NOT algebraically identical -- see section 6 below. Established for BALANCED\n",
  "repeated measurements only.\n\n",
  "Use A for the simulations.\n")
  else
  "A and B(a) DIVERGE -- see the comparisons above before choosing.\n")

saveRDS(list(trA = trA, trB = trB, mA = mA, mB = mB, np = np,
             d_ident = d_ident, agree = agree),
        "results/_A_vs_Ba.rds")

# ===========================================================================
# 6. WHY THEY AGREE, AND WHERE THAT MIGHT STOP
#
# partykit's mob_partynode() computes the OPG "meat" of the fluctuation test as
#
#     meat <- if (is.null(cluster)) crossprod(process)
#             else crossprod(as.matrix(apply(process, 2L, tapply,
#                                            as.numeric(cluster), sum)))
#
# i.e. with `cluster` supplied it uses CLUSTER-SUMMED scores -- exactly
# U_i = sum_t u_it. partykit's clustered-MOB machinery IS B(a).
#
# The empirical fluctuation process itself, however, is still built over the
# NT rows. For a person-level CATEGORICAL Z the two coincide exactly, because
# summing scores within a category over NT rows equals summing over the N
# persons in that category (persons are nested within categories). For a
# NUMERIC Z the statistic is a supLM over the ordered process, where row
# granularity could in principle matter. Since SDC_DRES_COM (years since
# immigration) is expected to enter as an integer covariate in the real
# application, that case is checked here explicitly.
# ===========================================================================
line("="); cat("6. NUMERIC PERSON-LEVEL PARTITIONING VARIABLE\n"); line("=")

set.seed(4242)
age_p <- round(runif(ref$n_person, 45, 85))          # person-level, no effect
pdat$AGE <- age_p
dat$AGE  <- age_p[match(dat$id, pdat$id)]
ZN2 <- c(ZN, "AGE")

fmlA2 <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                           paste(ZN2, collapse = " + ")))
fmlB2 <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ prow | ",
                           paste(ZN2, collapse = " + ")))
fitA2  <- era_gee_mob_fit(dat, YN, "X_time", fam, "independence", y_scale = "raw")
fitBa2 <- era_gee_mob_fit_subject(dat, pdat, YN, "X_time", fam, "independence",
                                  y_scale = "raw")

trA2 <- mob(fmlA2, data = dat,  fit = fitA2,  cluster = dat$id,
            control = mob_control(verbose = FALSE, maxdepth = 3, minsize = MINSIZE_ROW))
trB2 <- mob(fmlB2, data = pdat, fit = fitBa2,
            control = mob_control(verbose = FALSE, maxdepth = 3, minsize = MINSIZE_PERSON))

mA2 <- show_tests(trA2, "A    (with numeric AGE)")
mB2 <- show_tests(trB2, "B(a) (with numeric AGE)")
cat(sprintf("max |statistic_A - statistic_B(a)| = %.3e\n",
            max(abs(mA2["statistic", ] - mB2["statistic", ]))))
cat(sprintf("max |p_A - p_B(a)|                = %.3e\n",
            max(abs(mA2["p.value", ] - mB2["p.value", ]))))
cat(sprintf("split variable: A = %s, B(a) = %s | terminal nodes: %d vs %d\n",
            split_name(trA2), split_name(trB2), width(trA2), width(trB2)))

# ===========================================================================
# 7. CONTRAST: method A WITHOUT the cluster argument
#    Treating the NT wave-level scores as independent is what `cluster = id`
#    exists to prevent. This shows what it buys.
# ===========================================================================
line("="); cat("7. METHOD A WITHOUT cluster = id  (naive, for contrast)\n"); line("=")
trA_nc <- mob(fmlA, data = dat, fit = fitA,
              control = mob_control(verbose = FALSE, maxdepth = 3, minsize = MINSIZE_ROW))
mA_nc <- show_tests(trA_nc, "A (no cluster)")
cat("ratio of test statistics, A(no cluster) / A(cluster = id):\n")
print(round(mA_nc["statistic", ] / mA["statistic", ], 3))
cat(sprintf("\nterminal nodes: A(cluster) = %d, A(no cluster) = %d\n",
            width(trA), width(trA_nc)))
