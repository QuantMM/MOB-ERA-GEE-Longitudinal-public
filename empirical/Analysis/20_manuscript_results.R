# ---------------------------------------------------------------------------
# 20_manuscript_results.R
#
# Everything the manuscript Results section needs that is not already in the
# master drafts:
#
#   PART 1  Table 1 numbers: composition by immigrant status, including
#           education, income adequacy, housing tenure, and Wave-1 outcome
#           means, which were not previously tabulated by group.
#   PART 2  Figure 1: fitted four-domain trajectories by immigrant status at
#           the reference baseline-visit age.
#   PART 3  Figure 2: the full-sample slope-diagnostic MOB tree, drawn as a
#           schematic with node-specific per-wave changes.
#   PART 4  Figure 3: fitted trajectories for the prespecified racialization
#           strata within each immigration group. These strata are fixed in
#           advance and estimated identically in every group, so the figure
#           does not depend on tree shape or on group size.
#
# Figures are written as PDF (vector, for the submitted PDF) and TIFF at
# 300 dpi (for Frontiers figure upload). The values behind every figure are
# also written to CSV so that each plotted point can be checked.
#
# Run from the project root:
#     Rscript Analysis/20_manuscript_results.R
# ---------------------------------------------------------------------------

source("Analysis/era_gee_components.R")

OUT <- "DataPrep/out"
FIG <- "Draft/figures"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c(y_walk = "Walking", y_lsport = "Light activity",
         y_msport = "Moderate activity", y_ssport = "Strenuous activity")

bal <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
raw <- readRDS(file.path(OUT, "immigration_items_raw.rds"))

i  <- match(bal$entity_id, raw$entity_id)
fl <- raw$SDC_FIMM_COM[i]
yr <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[i]))
bal$imm <- factor(ifelse(fl == "2", "Non-immigrant",
                  ifelse(fl == "1", "Immigrant", NA)),
                  levels = c("Non-immigrant", "Immigrant"))
bal$grp3 <- factor(ifelse(fl == "2", "Non-immigrant",
                   ifelse(fl == "1" & !is.na(yr) & yr < 996,
                          ifelse(yr < 20, "Recent/mid-term", "Established"), NA)),
                   levels = c("Non-immigrant", "Recent/mid-term", "Established"))
bal <- bal[!is.na(bal$imm), ]
AGE_C <- 62.15

design <- function(d) { Ac <- d$z_age - AGE_C
  cbind(beta1 = 1, beta2 = d$X_time, gamma = Ac, delta = d$X_time * Ac) }
fit_one <- function(d) era_gee_components(d[, Y], design(d), d$id)

sup <- function(n) if (n > 0 && n < 6) "<6" else format(n, big.mark = ",", trim = TRUE)
pct <- function(n, d) sprintf("%.1f", 100 * n / d)

# ===========================================================================
cat("\n########## PART 1: Table 1 numbers ##########\n")
b0 <- bal[bal$wave == 0, ]
grps <- list(`Full sample` = rep(TRUE, nrow(b0)),
             `Non-immigrant` = b0$imm == "Non-immigrant",
             `Immigrant` = b0$imm == "Immigrant")

catrow <- function(lab, f) {
  cat(sprintf("  %-34s", lab))
  for (g in names(grps)) cat(sprintf(" %18s", f(b0[grps[[g]], ])))
  cat("\n")
}
cat(sprintf("  %-34s %18s %18s %18s\n", "", names(grps)[1], names(grps)[2], names(grps)[3]))
catrow("Participants, n", function(d) sup(nrow(d)))
catrow("Person-waves, n", function(d) sup(3 * nrow(d)))
catrow("Age at baseline visit, M (SD)",
       function(d) sprintf("%.2f (%.2f)", mean(d$z_age), sd(d$z_age)))
catrow("Age range", function(d) sprintf("%d-%d", min(d$z_age), max(d$z_age)))
catrow("Women, n (%)",
       function(d) sprintf("%s (%s)", sup(sum(d$z_sex == "Female")),
                           pct(sum(d$z_sex == "Female"), nrow(d))))
catrow("Racialized, n (%)",
       function(d) { k <- sum(d$z_racialized == "Racialized", na.rm = TRUE)
         sprintf("%s (%s)", sup(k), pct(k, sum(!is.na(d$z_racialized)))) })

catlevels <- function(var, label) {
  lv <- levels(b0[[var]])
  cat(sprintf("  %s\n", label))
  for (l in lv) catrow(paste0("   ", l),
    function(d) { k <- sum(d[[var]] == l, na.rm = TRUE)
      sprintf("%s (%s)", sup(k), pct(k, sum(!is.na(d[[var]])))) })
  catrow("   missing", function(d) sup(sum(is.na(d[[var]]))))
}
catlevels("z_educ4", "Education, n (%)")
catlevels("z_incneeds", "Income meets needs, n (%)")
catlevels("z_own", "Housing tenure, n (%)")

cat("\n  Wave 1 outcome means, M (SD)\n")
for (q in Y) catrow(paste0("   ", YL[[q]]),
  function(d) sprintf("%.3f (%.3f)", mean(d[[q]]), sd(d[[q]])))

cat("\n  Three-group version (for the supplement)\n")
for (g in levels(bal$grp3)) { s <- b0[!is.na(b0$grp3) & b0$grp3 == g, ]
  cat(sprintf("   %-16s n=%6s  age %.2f (%.2f)  %sF  %s racialized\n", g, sup(nrow(s)),
      mean(s$z_age), sd(s$z_age), pct(sum(s$z_sex == "Female"), nrow(s)),
      pct(sum(s$z_racialized == "Racialized", na.rm = TRUE), sum(!is.na(s$z_racialized))))) }

# ===========================================================================
cat("\n########## PART 2: Figure 1 ##########\n")
f_imm <- lapply(setNames(levels(bal$imm), levels(bal$imm)),
                function(g) fit_one(bal[bal$imm == g, ]))
fit_traj <- function(f) sapply(Y, function(q) f$B["beta1", q] + (0:2) * f$B["beta2", q])
traj1 <- lapply(f_imm, fit_traj)
for (g in names(traj1)) { cat("  ", g, "\n"); print(round(traj1[[g]], 4)) }

d1 <- do.call(rbind, lapply(names(traj1), function(g)
  data.frame(figure = 1, group = g, wave = 1:3, traj1[[g]], check.names = FALSE)))

COL_NON <- "#2a78d6"   # slot 1, blue: non-immigrants (both figures)
COL_IMM <- "#eb6834"   # slot 2, orange: immigrants (Figure 1)
draw_fig1 <- function() {
  op <- par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2.2, 0.8), mgp = c(2.3, 0.7, 0),
            cex.axis = 0.9, cex.lab = 0.95)
  on.exit(par(op))
  for (q in Y) {
    yv <- sapply(traj1, function(m) m[, q])
    rg <- range(yv); pad <- diff(rg) * 0.25 + 0.01
    plot(NA, xlim = c(0.9, 3.1), ylim = c(rg[1] - pad, rg[2] + pad),
         xaxt = "n", xlab = "CLSA assessment", ylab = "Model-implied scale points",
         main = YL[[q]], font.main = 1)
    axis(1, at = 1:3, labels = c("Baseline", "FU1", "FU2"))
    # Colors: categorical slots 1-2 of the validated dataviz reference palette
    # (blue, orange). Line type and marker shape are kept as secondary encoding.
    lines(1:3, yv[, "Non-immigrant"], lwd = 2, col = COL_NON)
    points(1:3, yv[, "Non-immigrant"], pch = 16, col = COL_NON, cex = 1.1)
    lines(1:3, yv[, "Immigrant"], lwd = 2, lty = 2, col = COL_IMM)
    points(1:3, yv[, "Immigrant"], pch = 17, col = COL_IMM, cex = 1.1)
    if (q == Y[1]) legend("bottomleft", c("Non-immigrant", "Immigrant"),
      lty = c(1, 2), pch = c(16, 17), col = c(COL_NON, COL_IMM),
      bty = "n", cex = 0.85, lwd = 2)
  }
}
pdf(file.path(FIG, "figure1.pdf"), width = 7.0, height = 6.0); draw_fig1(); dev.off()
tiff(file.path(FIG, "figure1.tiff"), width = 7.0, height = 6.0, units = "in",
     res = 300, compression = "lzw"); draw_fig1(); dev.off()
png(file.path(FIG, "figure1.png"), width = 7.0, height = 6.0, units = "in",
    res = 110); draw_fig1(); dev.off()

# ===========================================================================
cat("\n########## PART 3: Figure 2 ##########\n")
mob <- readRDS(file.path(OUT, "_mob_final_model.rds"))
sl  <- mob[["D1 slope"]]
node <- sl$node; off <- sl$offset

# The saved tree was fitted on the balanced file WITHOUT dropping participants
# whose Immigration Flag is missing, so node membership is indexed against that
# row ordering. Reload it unfiltered; using the filtered object here would shift
# every index and silently mis-assign nodes.
balF <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
iF   <- match(balF$entity_id, raw$entity_id)
flF  <- raw$SDC_FIMM_COM[iF]
balF$imm <- factor(ifelse(flF == "2", "Non-immigrant",
                   ifelse(flF == "1", "Immigrant", NA)),
                   levels = c("Non-immigrant", "Immigrant"))
stopifnot(length(node) == nrow(balF), nrow(off) == nrow(balF))
bal2 <- balF
nodes <- sort(unique(node[!is.na(node)]))
ninfo <- lapply(nodes, function(nd) {
  ix <- which(!is.na(node) & node == nd); s <- bal2[ix, ]
  yy <- as.matrix(s[, Y]) - off[ix, , drop = FALSE]
  m  <- era_gee_components(yy, cbind(beta1 = 1, beta2 = s$X_time), s$id)
  list(n = length(unique(s$id)), b2 = m$B["beta2", Y],
       pF = 100 * mean(s$z_sex == "Female", na.rm = TRUE),
       pI = 100 * mean(s$imm == "Immigrant", na.rm = TRUE))
})
names(ninfo) <- as.character(nodes)
for (k in names(ninfo)) cat(sprintf("  node %s: n=%s  %.0f%%F  %.0f%%imm  b2 %s\n", k,
  sup(ninfo[[k]]$n), ninfo[[k]]$pF, ninfo[[k]]$pI,
  paste(sprintf("%+.4f", ninfo[[k]]$b2), collapse = " ")))

draw_fig2 <- function() {
  op <- par(mar = c(0.5, 0.5, 0.5, 0.5)); on.exit(par(op))
  plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "")
  oval <- function(x, y, lab) {
    rect(x - 14, y - 4.5, x + 14, y + 4.5, col = "grey92", border = "grey30")
    text(x, y, lab, cex = 0.85) }
  nd <- function(x, y, title, info) {
    rect(x - 16, y - 14, x + 16, y + 14, col = "white", border = "grey30")
    text(x, y + 10, title, cex = 0.78, font = 2)
    text(x, y - 2, info, cex = 0.72) }
  fmt <- function(k) paste(sprintf("%s %+.3f", c("Walk", "Light", "Mod.", "Stren."),
                                   ninfo[[k]]$b2), collapse = "\n")
  # Terminal boxes are placed on disjoint horizontal bands so that the female
  # terminal node cannot collide with the male sub-tree.
  oval(50, 93, "Sex")
  segments(50, 88.5, 16, 74); segments(50, 88.5, 70, 79.5)
  text(30, 83, "Female", cex = 0.75); text(63, 86, "Male", cex = 0.75)
  nd(16, 60, sprintf("Women (n = %s)", sup(ninfo[["2"]]$n)), fmt("2"))
  oval(70, 75, "Immigrant status")
  segments(70, 70.5, 46, 42); segments(70, 70.5, 84, 42)
  text(52, 58, "Non-immigrant", cex = 0.72); text(82, 58, "Immigrant", cex = 0.72)
  nd(46, 28, sprintf("Non-immigrant men (n = %s)", sup(ninfo[["4"]]$n)), fmt("4"))
  nd(84, 28, sprintf("Immigrant men (n = %s)", sup(ninfo[["5"]]$n)), fmt("5"))
  text(50, 5, paste("Entries are node-specific change per data-collection wave",
                    "in scale points,\nat the reference baseline-visit age of 62.15 years."),
       cex = 0.72)
}
pdf(file.path(FIG, "figure2.pdf"), width = 7.0, height = 5.2); draw_fig2(); dev.off()
tiff(file.path(FIG, "figure2.tiff"), width = 7.0, height = 5.2, units = "in",
     res = 300, compression = "lzw"); draw_fig2(); dev.off()
png(file.path(FIG, "figure2.png"), width = 7.0, height = 5.2, units = "in",
    res = 110); draw_fig2(); dev.off()

# ===========================================================================
cat("\n########## PART 4: Figure 3 ##########\n")
cells <- expand.grid(grp = levels(bal$grp3), rac = c("White", "Racialized"),
                     stringsAsFactors = FALSE)
f_cell <- list(); traj3 <- list()
for (k in seq_len(nrow(cells))) {
  g <- cells$grp[k]; r <- cells$rac[k]
  s <- bal[!is.na(bal$grp3) & bal$grp3 == g & !is.na(bal$z_racialized) & bal$z_racialized == r, ]
  np <- length(unique(s$id))
  nm <- paste(g, r, sep = " / ")
  if (np < 30) { cat(sprintf("  %-34s n=%s  skipped\n", nm, sup(np))); next }
  f <- fit_one(s); f_cell[[nm]] <- f; traj3[[nm]] <- fit_traj(f)
  cat(sprintf("  %-34s n=%6s  b1 %s | b2 %s\n", nm, sup(np),
              paste(sprintf("%.3f", f$B["beta1", Y]), collapse = " "),
              paste(sprintf("%+.3f", f$B["beta2", Y]), collapse = " ")))
}
d3 <- do.call(rbind, lapply(names(traj3), function(nm)
  data.frame(figure = 3, group = nm, wave = 1:3, traj3[[nm]], check.names = FALSE)))
write.csv(rbind(d1, d3), file.path(FIG, "figure_data.csv"), row.names = FALSE)

# Categorical slots 1-3 of the validated dataviz reference palette, which pass the
# all-pairs checks. Legend order follows slot order; established immigrants (most
# immigrants) take slot 2 so the immigrant color carries over from Figure 1.
# Slot 3 (aqua) is below 3:1 contrast on white, so line type, marker fill, and the
# caption carry the encoding as well.
COLG <- c(`Non-immigrant` = "#2a78d6", `Established` = "#eb6834", `Recent/mid-term` = "#1baf7a")
draw_fig3 <- function() {
  op <- par(mfrow = c(2, 2), mar = c(3.6, 3.8, 2.2, 0.8), mgp = c(2.3, 0.7, 0),
            cex.axis = 0.9, cex.lab = 0.95); on.exit(par(op))
  for (q in Y) {
    yv <- sapply(traj3, function(m) m[, q])
    rg <- range(yv); pad <- diff(rg) * 0.28 + 0.01
    plot(NA, xlim = c(0.9, 3.1), ylim = c(rg[1] - pad, rg[2] + pad), xaxt = "n",
         xlab = "CLSA assessment", ylab = "Model-implied scale points",
         main = YL[[q]], font.main = 1)
    axis(1, at = 1:3, labels = c("Baseline", "FU1", "FU2"))
    for (nm in names(traj3)) {
      g <- sub(" / .*$", "", nm); r <- sub("^.* / ", "", nm)
      lines(1:3, traj3[[nm]][, q], col = COLG[[g]], lwd = 2,
            lty = if (r == "White") 1 else 2)
      points(1:3, traj3[[nm]][, q], col = COLG[[g]],
             pch = if (r == "White") 16 else 1)
    }
    # the two legends are placed in different panels so that neither sits on
    # top of a fitted line
    if (q == Y[1]) legend("topleft", c("White", "Racialized"), lty = c(1, 2),
                          pch = c(16, 1), bty = "n", cex = 0.72)
    if (q == Y[4]) legend("bottomleft", names(COLG), col = COLG, lwd = 2,
                          bty = "n", cex = 0.72)
  }
}
pdf(file.path(FIG, "figure3.pdf"), width = 7.0, height = 6.0); draw_fig3(); dev.off()
tiff(file.path(FIG, "figure3.tiff"), width = 7.0, height = 6.0, units = "in",
     res = 300, compression = "lzw"); draw_fig3(); dev.off()
png(file.path(FIG, "figure3.png"), width = 7.0, height = 6.0, units = "in",
    res = 110); draw_fig3(); dev.off()

cat("\nwritten to", FIG, ":\n  figure1/2/3 .pdf and .tiff, figure_data.csv\n")
