# ---------------------------------------------------------------------------
# 01_synthetic_example.R
#
# One complete, substantively interpretable end-to-end example of the agreed
# MOB + ERA-GEE base model, on a synthetic dataset built to resemble the
# intended CLSA application rather than an abstract simulation.
#
# DESIGN OF THE TRUTH (deliberately simple and obvious)
#   * baseline levels are IDENTICAL across immigrant groups
#   * recent immigrants decline substantially faster on ALL FOUR outcomes
#   * no other partitioning variable generates any heterogeneity
#   -> MOB should recover immigration status, separating "recent" from the rest,
#      and it should do so because of SLOPE instability, not intercept instability.
#
# Run from the project root:
#     Rscript run/01_worked_example.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("R/era_gee_fixed.R")
source("R/synthetic_data.R")

# ===========================================================================
# 1. THE REFERENCE CASE  (generator lives in R/synthetic_data.R so that
#    every script in this folder works on byte-identical data)
# ===========================================================================
ref      <- make_reference_data()
dat      <- ref$dat
person   <- ref$person
YN       <- ref$YN
ZN       <- ref$ZN
Q        <- ref$Q
N_PERSON <- ref$n_person
N_WAVE   <- ref$n_wave

TRUE_B1        <- ref$truth$b1
TRUE_B2_OTHERS <- ref$truth$b2_others
TRUE_B2_RECENT <- ref$truth$b2_recent

# ===========================================================================
# 4. REPORT: data structure and group sizes
# ===========================================================================
line <- function(ch = "-") cat(strrep(ch, 78), "\n")
line("="); cat("1. GENERATED DATA STRUCTURE\n"); line("=")
cat(sprintf("persons = %d | waves = %d (t = 0,1,2) | rows = %d | outcomes = %d\n\n",
            N_PERSON, N_WAVE, nrow(dat), Q))
cat("First rows (person-wave layout):\n")
print(head(dat[, c("id","wave","X_time",YN,"IMM","SEX")], 6), digits = 3)

cat("\nPerson-level group sizes (partitioning variables):\n")
for (z in ZN) {
  tb <- table(person[[z]])
  cat(sprintf("  %-8s %s\n", z,
      paste(sprintf("%s=%d (%.1f%%)", names(tb), tb, 100*tb/N_PERSON), collapse = "  ")))
}

cat("\nObserved outcome means by wave and immigrant group:\n")
agg <- aggregate(dat[, YN], by = list(IMM = dat$IMM, wave = dat$wave), mean)
print(agg[order(agg$IMM, agg$wave), ], row.names = FALSE, digits = 3)

# ===========================================================================
# 5. REPORT: the truth used in generation
# ===========================================================================
line("="); cat("2. TRUE PARAMETERS USED IN GENERATION\n"); line("=")
truth <- rbind(`intercept (all groups)` = TRUE_B1,
               `slope: non / established` = TRUE_B2_OTHERS,
               `slope: recent`            = TRUE_B2_RECENT)
print(round(truth, 3))
cat("\nBaselines are IDENTICAL across immigrant groups by construction.\n")
cat("Only the slopes differ, and only for 'recent'.\n")
cat("SEX, ETHN, EDU, INCNEED, HOMEOWN, URBAN have NO effect on either.\n")

# ===========================================================================
# 6. ROOT MODEL
# ===========================================================================
line("="); cat("3. FITTED ROOT MODEL (all persons pooled)\n"); line("=")
t0 <- Sys.time()
root <- ERA_GEE_fixed(dat[, YN], dat$X_time, rep("gaussian", Q),
                      corstr = "independence", y_scale = "raw", verbose = TRUE)
cat(sprintf("elapsed %.1fs | B is %d x %d | estfun is %d x %d | colSums(estfun) max |.| = %.2e\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            nrow(root$B), ncol(root$B), nrow(root$estfun), ncol(root$estfun),
            max(abs(colSums(root$estfun)))))
cat("\nRoot coefficient matrix B (2 x Q):\n"); print(round(root$B, 4))
cat("\nThe root slopes are a mixture: 12% of persons decline steeply, 88% mildly.\n")

# ===========================================================================
# 7. MOB TREE
# ===========================================================================
line("="); cat("4. MOB-ERA-GEE TREE\n"); line("=")
fitf <- era_gee_mob_fit(dat, YN, "X_time", rep("gaussian", Q),
                        corstr = "independence", y_scale = "raw")
fml <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                         paste(ZN, collapse = " + ")))
t0 <- Sys.time()
tr <- mob(fml, data = dat, fit = fitf, cluster = dat$id,
          control = mob_control(verbose = TRUE, maxdepth = 3, minsize = 450))
cat(sprintf("\ntree fitted in %.1f min | terminal nodes = %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), width(tr)))
cat("\n"); print(tr)

# ===========================================================================
# 8. NODE ESTIMATES VS TRUTH
# ===========================================================================
line("="); cat("5. TERMINAL-NODE ESTIMATES vs TRUTH\n"); line("=")
leaves <- nodeids(tr, terminal = TRUE)
memb   <- predict(tr, type = "node")

cmp <- do.call(rbind, lapply(leaves, function(nd) {
  d  <- dat[memb == nd, ]
  m  <- ERA_GEE_fixed(d[, YN], d$X_time, rep("gaussian", Q),
                      corstr = "independence", y_scale = "raw")
  gp <- sort(unique(as.character(d$IMM[!duplicated(d$id)])))
  data.frame(node = nd,
             n_person = length(unique(d$id)), n_row = nrow(d),
             IMM_groups = paste(gp, collapse = " + "),
             outcome = YN,
             est_intercept = round(m$B["intercept", ], 3),
             true_intercept = round(TRUE_B1, 3),
             est_slope = round(m$B["slope", ], 3),
             true_slope = round(if (identical(gp, "recent")) TRUE_B2_RECENT
                                else TRUE_B2_OTHERS, 3),
             row.names = NULL)
}))
print(cmp, row.names = FALSE)

# ===========================================================================
# 9. TRAJECTORY PLOT
# ===========================================================================
png("results/fig_synthetic_trajectories.png", width = 1100, height = 850, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4.2, 3, 1), oma = c(7, 0, 3, 0))
cols <- c("#1b6ca8", "#d1495b", "#2e8b57", "#8e6c8a")
for (q in seq_len(Q)) {
  yl <- range(sapply(leaves, function(nd)
        tapply(dat[memb == nd, YN[q]], dat$wave[memb == nd], mean)))
  yl <- yl + c(-.12, .12) * diff(yl)
  plot(NA, xlim = c(-0.1, 2.1), ylim = yl, xaxt = "n",
       xlab = "wave (t)", ylab = YN[q], main = YN[q])
  axis(1, at = 0:2)
  for (k in seq_along(leaves)) {
    nd <- leaves[k]; d <- dat[memb == nd, ]
    mo <- tapply(d[[YN[q]]], d$wave, mean)
    b  <- cmp[cmp$node == nd & cmp$outcome == YN[q], ]
    points(0:2, mo, pch = 19, col = cols[k], cex = 1.15)
    abline(a = b$est_intercept, b = b$est_slope, col = cols[k], lwd = 2)
    abline(a = TRUE_B1[q], b = b$true_slope, col = cols[k], lwd = 1, lty = 3)
  }
}
labs <- sapply(leaves, function(nd) {
  d  <- dat[memb == nd, ]
  gp <- sort(unique(as.character(d$IMM[!duplicated(d$id)])))
  sprintf("node %d: IMM = %s   (n = %d persons)", nd,
          paste(gp, collapse = " + "), length(unique(d$id)))
})
nl <- length(leaves)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", inset = 0.005, legend = labs, col = cols[seq_len(nl)],
       lwd = 2, pch = 19, bty = "n", cex = .9, ncol = 1,
       title = "terminal nodes", title.adj = 0)
legend("bottomright", inset = c(0.04, 0.005),
       legend = c("points = observed wave means", "solid  = fitted trajectory",
                  "dotted = true trajectory"),
       bty = "n", cex = .85)
mtext("MOB-ERA-GEE on synthetic data: fitted vs true trajectories",
      outer = TRUE, line = -2, cex = 1.1, font = 2)
par(op); invisible(dev.off())
cat("\nplot written to results/fig_synthetic_trajectories.png\n")

saveRDS(list(dat = dat, tree = tr, root = root, cmp = cmp),
        "results/_synthetic_example.rds")
