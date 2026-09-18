# ---------------------------------------------------------------------------
# 07_attrition_and_selection.R
#
# Two descriptive checks, both requested before settling the analytic sample.
#
# (A) ATTRITION. Compare BASELINE outcome distributions between participants
#     with complete and incomplete follow-up. Comparing baseline values is the
#     only clean comparison, since people missing a follow-up wave have no
#     outcome there by definition.
#
#     No significance tests are reported. At N ~ 30,000 a trivial difference is
#     "significant", so p-values would mislead rather than inform. This does NOT
#     establish MCAR. The claim it supports is narrower: participants with
#     complete and incomplete follow-up showed similar baseline distributions on
#     the study outcomes.
#
# (B) SELECTION INTO THE BALANCED SAMPLE. Three-wave completeness cross-tabulated
#     with immigration status, racialization, sex and baseline age, to check in
#     particular whether more recent immigrants are disproportionately lost when
#     the analysis is restricted to complete cases.
#
# Run from the project root:
#     Rscript DataPrep/07_attrition_and_selection.R
# Writes DataPrep/out/attrition_selection.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

OUT <- "DataPrep/out"
a   <- readRDS(file.path(OUT, "analytic_long_all.rds"))
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c("Walking", "Light exercise", "Moderate exercise", "Strenuous exercise")

con <- file(file.path(OUT, "attrition_selection.md"), open = "wt")
w   <- function(...) writeLines(paste0(...), con)
supp <- function(x, thr = 6) {
  x <- as.integer(x); o <- format(x, big.mark = ",", trim = TRUE)
  s <- !is.na(x) & x > 0 & x < thr; o[s] <- paste0("<", thr)
  if (sum(s) == 1L && length(x) > 1L) {
    c2 <- which(!s & !is.na(x)); if (length(c2)) o[c2[which.min(x[c2])]] <- paste0("<", thr)
  }
  o
}
mdtab <- function(m, rn = "") {
  w("| ", rn, " | ", paste(colnames(m), collapse = " | "), " |")
  w("|", paste(rep("---", ncol(m) + 1), collapse = "|"), "|")
  for (i in seq_len(nrow(m))) w("| ", rownames(m)[i], " | ", paste(m[i, ], collapse = " | "), " |")
  w("")
}

# person-level frame, baseline row
b <- a[a$wave == 0, ]
b$complete <- factor(ifelse(b$complete3, "Complete (3 waves)", "Incomplete (<3 waves)"),
                     levels = c("Complete (3 waves)", "Incomplete (<3 waves)"))

w("# Attrition and selection into the balanced sample")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `DataPrep/07_attrition_and_selection.R`")
w("")
w("Counts below 6 are suppressed per the CLSA Publication Policy.")
w("")
w("---")
w("")

# ===========================================================================
w("## A. Attrition — baseline outcomes by follow-up completeness")
w("")
w("Participants observed at baseline with all four outcomes: **",
  format(nrow(b), big.mark = ","), "**.")
w("")
tb <- table(b$complete)
m <- cbind(persons = supp(as.integer(tb)),
           `%` = sprintf("%.1f", 100 * as.integer(tb) / nrow(b)))
rownames(m) <- names(tb)
mdtab(m, "follow-up")
w("Nonparticipation is predominantly monotonic: of the ",
  format(nrow(b), big.mark = ","), " baseline participants, only ",
  supp(sum(!is.na(b$entity_id) & b$n_waves == 2 &
           sapply(b$entity_id, function(i) identical(sort(a$wave[a$entity_id == i]), c(0L, 2L))))),
  " were absent at Wave 1 but returned at Wave 2.")
w("")
w("### Baseline outcome distributions")
w("")
m <- t(sapply(seq_along(Y), function(i) {
  s <- split(b[[Y[i]]], b$complete)
  c(sprintf("%.2f (%.2f)", mean(s[[1]], na.rm=TRUE), sd(s[[1]], na.rm=TRUE)),
    sprintf("%.2f (%.2f)", mean(s[[2]], na.rm=TRUE), sd(s[[2]], na.rm=TRUE)),
    sprintf("%+.2f", mean(s[[2]], na.rm=TRUE) - mean(s[[1]], na.rm=TRUE)),
    sprintf("%+.3f", (mean(s[[2]], na.rm=TRUE) - mean(s[[1]], na.rm=TRUE)) /
                      sd(b[[Y[i]]], na.rm = TRUE)))
}))
colnames(m) <- c("Complete: M (SD)", "Incomplete: M (SD)", "difference", "std. difference")
rownames(m) <- YL
mdtab(m, "baseline outcome")
smd <- sapply(Y, function(v) {
  s <- split(b[[v]], b$complete)
  (mean(s[[2]], na.rm=TRUE) - mean(s[[1]], na.rm=TRUE)) / sd(b[[v]], na.rm = TRUE)
})
w("No significance tests are reported: at N = ", format(nrow(b), big.mark = ","),
  " they would flag differences of no")
w("practical size. Standardised differences are used instead.")
w("")
w("**The differences are small but not negligible, and they are systematic.** All")
w("four point the same way: participants who did not complete follow-up were")
w("*less* physically active at baseline. The gradient tracks intensity — it is")
w("smallest for light exercise (", sprintf("%.3f", smd[2]), ") and largest for")
w("strenuous exercise (", sprintf("%.3f", smd[4]), "), with walking (",
  sprintf("%.3f", smd[1]), ") and moderate")
w("exercise (", sprintf("%.3f", smd[3]), ") in between. Three of the four exceed the 0.1")
w("threshold often used for negligible imbalance, though all are well below 0.25.")
w("")
w("This is the pattern health-related attrition would produce, and it is")
w("consistent with the age gradient in Section B. It should be described")
w("accurately rather than as \"highly similar\": the honest statement is that")
w("baseline outcome differences between completers and non-completers were small")
w("(standardised differences ", sprintf("%.2f", min(abs(smd))), " to ",
  sprintf("%.2f", max(abs(smd))), "), with non-completers somewhat less active,")
w("most visibly for strenuous activity.")
w("")
w("### Response distribution, for completeness")
w("")
for (i in seq_along(Y)) {
  tt <- table(b[[Y[i]]], b$complete)
  m <- apply(tt, 2, function(x) sprintf("%s (%.1f%%)", supp(x), 100 * x / sum(x)))
  rownames(m) <- c("1 Never", "2 Seldom", "3 Sometimes", "4 Often")
  w("**", YL[i], "**"); w(""); mdtab(m, "baseline response")
}

# ===========================================================================
w("## B. Who is lost when the sample is restricted to three complete waves")
w("")
w("Percentage of each baseline group retained in the balanced sample. The")
w("question this is meant to answer is whether the groups the study is about are")
w("disproportionately lost by a complete-case restriction.")
w("")
retain_by <- function(v, lab) {
  if (all(is.na(b[[v]]))) return(invisible(NULL))
  tt <- table(b[[v]], b$complete)
  m <- cbind(complete = supp(tt[, 1]), incomplete = supp(tt[, 2]),
             `% retained` = sprintf("%.1f", 100 * tt[, 1] / rowSums(tt)),
             `n baseline` = supp(rowSums(tt)))
  rownames(m) <- rownames(tt)
  w("**", lab, "**"); w(""); mdtab(m, lab)
}
retain_by("z_imm_cut20", "Immigration status (20-year threshold)")
retain_by("z_immigrant", "Immigrant vs non-immigrant")
retain_by("z_racialized", "Racialization")
retain_by("z_g4_cut20", "Racialization x immigration (20-year threshold)")
retain_by("z_sex", "Sex")
retain_by("z_educ4", "Education")

w("**Baseline age**")
w("")
b$agegrp <- cut(b$z_age, breaks = c(44, 54, 64, 74, 100),
                labels = c("45-54", "55-64", "65-74", "75+"))
tt <- table(b$agegrp, b$complete)
m <- cbind(complete = supp(tt[, 1]), incomplete = supp(tt[, 2]),
           `% retained` = sprintf("%.1f", 100 * tt[, 1] / rowSums(tt)),
           `n baseline` = supp(rowSums(tt)))
rownames(m) <- rownames(tt)
mdtab(m, "baseline age")
m <- rbind(c(sprintf("%.1f (%.1f)", mean(b$z_age[b$complete3], na.rm=TRUE),
                     sd(b$z_age[b$complete3], na.rm=TRUE)),
             sprintf("%.1f (%.1f)", mean(b$z_age[!b$complete3], na.rm=TRUE),
                     sd(b$z_age[!b$complete3], na.rm=TRUE))))
colnames(m) <- c("Complete: M (SD)", "Incomplete: M (SD)"); rownames(m) <- "baseline age"
mdtab(m, "")

w("**Years since immigration, immigrants only**")
w("")
im <- b[!is.na(b$z_yrs_imm), ]
m <- rbind(c(sprintf("%.1f (%.1f)", mean(im$z_yrs_imm[im$complete3]),
                     sd(im$z_yrs_imm[im$complete3])),
             sprintf("%.1f (%.1f)", mean(im$z_yrs_imm[!im$complete3]),
                     sd(im$z_yrs_imm[!im$complete3])),
             sprintf("%.1f", 100 * sum(im$complete3) / nrow(im))))
colnames(m) <- c("Complete: M (SD)", "Incomplete: M (SD)", "% retained")
rownames(m) <- "years since immigration"
mdtab(m, "")

# --- summary of the selection gradient ------------------------------------
r <- function(v, lv) { tt <- table(b[[v]], b$complete); 100 * tt[lv, 1] / sum(tt[lv, ]) }
w("### Reading Section B")
w("")
w("Retention is not uniform, and the gradient runs against the groups of")
w("substantive interest, though modestly:")
w("")
w("- immigration status: non-immigrant ", sprintf("%.1f%%", r("z_imm_cut20","non-immigrant")),
  " · established ", sprintf("%.1f%%", r("z_imm_cut20","established")),
  " · recent ", sprintf("%.1f%%", r("z_imm_cut20","recent")))
w("- racialization: White ", sprintf("%.1f%%", r("z_racialized","White")),
  " · Racialized ", sprintf("%.1f%%", r("z_racialized","Racialized")))
w("- the combined variable is the sharpest: **Racialized + recent immigrant is")
w("  retained at ", sprintf("%.1f%%", r("z_g4_cut20","3 Racialized + recent immigrant")),
  "**, the lowest of the four categories, against ",
  sprintf("%.1f%%", r("z_g4_cut20","1 White")), " for White")
ta <- table(b$agegrp, b$complete)
ra <- 100 * ta[, 1] / rowSums(ta)
w("- baseline age dominates everything else: 45-54 ", sprintf("%.1f%%", ra[["45-54"]]),
  " versus 75+ ", sprintf("%.1f%%", ra[["75+"]]))
w("- education runs the same way: < secondary ",
  sprintf("%.1f%%", r("z_educ4","< secondary")), " versus post-secondary degree ",
  sprintf("%.1f%%", r("z_educ4","Post-secondary degree")))
w("")
w("Two things follow. First, the immigrant and racialization gradients are")
w("partly compositional: both groups are correlated with age and education, which")
w("show much stronger retention gradients. Second, the differential is real but")
w("small in absolute terms — the balanced sample still contains ",
  supp(sum(b$complete3 & b$z_imm_cut20 == "recent", na.rm = TRUE)),
  " more recent")
w("immigrants and ",
  supp(sum(b$complete3 & b$z_g4_cut20 == "3 Racialized + recent immigrant", na.rm = TRUE)),
  " racialized more recent immigrants, both well above any")
w("threshold at which node-level estimation would be a concern.")
w("")
w("A complete-case primary analysis is defensible on these numbers, provided the")
w("selection is stated rather than passed over: participants lost to follow-up")
w("were older, less educated, less physically active at baseline, and somewhat")
w("more likely to be racialized or more recent immigrants.")
w("")

# ===========================================================================
w("## C. Group sizes in the balanced sample")
w("")
w("What is actually available for the primary analysis if it is restricted to")
w("participants complete at all three waves.")
w("")
bb <- b[b$complete3, ]
for (v in c("z_imm_cut20", "z_racialized", "z_g4_cut20", "z_sex")) {
  tb <- table(bb[[v]])
  m <- cbind(n = supp(as.integer(tb)), `%` = sprintf("%.1f", 100 * as.integer(tb) / nrow(bb)))
  rownames(m) <- names(tb)
  w("**", v, "**"); w(""); mdtab(m, "")
}
tt <- table(bb$z_g4_cut20, bb$z_sex)
m <- apply(tt, 2, supp); rownames(m) <- rownames(tt)
w("**Racialization x immigration x sex** (smallest cell governs what MOB can find)")
w(""); mdtab(m, "")

close(con)
cat("written:", file.path(OUT, "attrition_selection.md"), "\n")
