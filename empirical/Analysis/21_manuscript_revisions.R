# ---------------------------------------------------------------------------
# 21_manuscript_revisions.R
#
# Numbers needed to revise the Results section (14 Sep 2026 review):
#
#   PART 1  Foundational model refitted on the Flag-determinate balanced sample,
#           so every Results number shares one denominator (the group analyses
#           already exclude participants whose Immigration Flag is missing;
#           reporting both totals would disclose a count below 6).
#   PART 2  Flag-determinate size of the all-available person-wave sample.
#   PART 3  Racialized-minus-White contrast compared ACROSS groups:
#           D = (Rac - White)_immigrant group - (Rac - White)_non-immigrant.
#           Groups and strata are disjoint, so Var(D) is the sum of four
#           stratum variances. Block Wald tests on 4 df.
#   PART 4  Instability tests at every inner node of the full-sample slope tree.
#
# The reference age is the same constant used in Analysis/15-20.
# Counts below 6 are suppressed. Run from the project root:
#     Rscript Analysis/21_manuscript_revisions.R
# Writes (local only): DataPrep/out/MANUSCRIPT_REVISIONS.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(partykit))
source("Analysis/era_gee_components.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c(y_walk = "Walking", y_lsport = "Light", y_msport = "Moderate", y_ssport = "Strenuous")
CMP <- c("beta1", "beta2", "gamma", "delta")

balF <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
alls <- readRDS(file.path(OUT, "analytic_long_all.rds"))
raw  <- readRDS(file.path(OUT, "immigration_items_raw.rds"))
AGE_C <- mean(balF$z_age[balF$wave == 0])          # identical to Analysis/15-20

sup <- function(n) if (n > 0 && n < 6) "<6" else format(n, big.mark = ",", trim = TRUE)
fp  <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))

add_grp <- function(d) {
  i  <- match(d$entity_id, raw$entity_id)
  fl <- raw$SDC_FIMM_COM[i]; yr <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[i]))
  d$imm  <- ifelse(fl == "2", "Non-immigrant", ifelse(fl == "1", "Immigrant", NA))
  d$grp3 <- ifelse(fl == "2", "Non-immigrant",
            ifelse(fl == "1" & !is.na(yr) & yr < 996,
                   ifelse(yr < 20, "Recent/mid-term", "Established"), NA))
  d
}
bal  <- add_grp(balF); bal <- bal[!is.na(bal$imm), ]
alls <- add_grp(alls); allsD <- alls[!is.na(alls$imm), ]

design  <- function(d, Ac = d$z_age - AGE_C) cbind(beta1 = 1, beta2 = d$X_time, gamma = Ac, delta = d$X_time * Ac)
fit_one <- function(d) era_gee_components(d[, Y], design(d), d$id)

con <- file(file.path(OUT, "MANUSCRIPT_REVISIONS.md"), open = "wt")
w <- function(...) { s <- paste0(...); writeLines(s, con); cat(s, "\n") }
w("# Numbers for the Results revision"); w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " - `Analysis/21_manuscript_revisions.R`"); w("")
w("Reference age: ", sprintf("%.4f", AGE_C)); w("")

# ===========================================================================
w("## PART 1: foundational model, Flag-determinate balanced sample"); w("")
b0 <- bal[bal$wave == 0, ]
w("- persons ", sup(nrow(b0)), ", person-waves ", sup(nrow(bal)),
  "; mean (SD) baseline-visit age ", sprintf("%.2f (%.2f)", mean(b0$z_age), sd(b0$z_age)))
cr <- cor(b0[, Y]); w("- Wave-1 correlations: range ", sprintf("%.3f to %.3f", min(cr[upper.tri(cr)]), max(cr[upper.tri(cr)])))
w("- Wave-1 means: ", paste(sprintf("%s %.3f", YL[Y], colMeans(b0[, Y])), collapse = ", "))
w("- Observed means W1/W2/W3: ", paste(sapply(Y, function(q) paste0(YL[[q]], " ",
  paste(sprintf("%.3f", sapply(0:2, function(k) mean(bal[[q]][bal$wave == k]))), collapse = "/"))), collapse = "; "))

f <- fit_one(bal)
w(""); w("| component | ", paste(YL[Y], collapse = " | "), " |"); w("|---|---|---|---|---|")
for (k in CMP) w("| ", k, " | ", paste(sprintf("%+.5f (%.5f)", f$B[k, ], f$SE[k, ]), collapse = " | "), " |")
w("")
for (k in c("beta2", "gamma", "delta")) { t <- era_wald(f, paste0(k, ":", Y))
  w("- joint ", k, " = 0: W = ", sprintf("%.2f", t["W"]), ", p = ", format(t["p"], digits = 3)) }
w("- single-outcome p: ", paste(sapply(CMP, function(k) paste0(k, " ",
  paste(fp(2 * pnorm(-abs(f$B[k, ] / f$SE[k, ]))), collapse = "/"))), collapse = "; "))
for (off in c(-10, 10)) w("- per-wave change at age ", sprintf("%.2f", AGE_C + off), ": ",
  paste(sprintf("%+.3f", f$B["beta2", ] + off * f$B["delta", ]), collapse = ", "))
w("- 20*delta (72.15 minus 52.15): ", paste(sprintf("%+.3f (%.3f)", 20 * f$B["delta", ], 20 * f$SE["delta", ]), collapse = ", "))

# equal-change check: saturated wave means, contrast mu0 - 2 mu1 + mu2 = -2 b_w2 + b_w3
Fs <- cbind(w1 = 1, w2 = as.numeric(bal$wave == 1), w3 = as.numeric(bal$wave == 2))
fs <- era_gee_components(bal[, Y], Fs, bal$id)
L  <- matrix(0, 4, length(fs$theta)); colnames(L) <- rownames(fs$V)
for (j in seq_along(Y)) { L[j, paste0("w2:", Y[j])] <- -2; L[j, paste0("w3:", Y[j])] <- 1 }
dd <- L %*% fs$theta; Vd <- L %*% fs$V %*% t(L)
Weq <- as.numeric(t(dd) %*% solve(Vd) %*% dd)
w("- equal-change joint test: W = ", sprintf("%.2f", Weq), ", p = ", format(pchisq(Weq, 4, lower.tail = FALSE), digits = 3))
w("- difference of changes by outcome: ", paste(sprintf("%s %+.3f (p %s)", YL[Y], dd, fp(2 * pnorm(-abs(dd / sqrt(diag(Vd)))))), collapse = ", "))

# ===========================================================================
w(""); w("## PART 2: all-available person-wave sample, Flag-determinate"); w("")
w("- persons ", sup(length(unique(allsD$id))), ", person-waves ", sup(nrow(allsD)))

# ===========================================================================
w(""); w("## PART 3: racialized-minus-White contrast, compared across groups"); w("")
strat <- function(g, r) bal[!is.na(bal$grp3) & bal$grp3 == g & !is.na(bal$z_racialized) & bal$z_racialized == r, ]
G3 <- c("Non-immigrant", "Recent/mid-term", "Established")
ct <- lapply(setNames(G3, G3), function(g) { fw <- fit_one(strat(g, "White")); fr <- fit_one(strat(g, "Racialized"))
  list(d = fr$theta - fw$theta, V = fr$V + fw$V, nm = rownames(fw$V)) })
# immigrants pooled (both duration groups), for a single immigrant-vs-non summary
stI <- function(r) bal[bal$imm == "Immigrant" & !is.na(bal$z_racialized) & bal$z_racialized == r, ]
fwI <- fit_one(stI("White")); frI <- fit_one(stI("Racialized"))
ct[["Immigrant (all)"]] <- list(d = frI$theta - fwI$theta, V = frI$V + fwI$V, nm = rownames(fwI$V))

nm <- ct[[1]]$nm
wald <- function(d, V, j) { W <- as.numeric(t(d[j]) %*% solve(V[j, j]) %*% d[j]); c(W = W, p = pchisq(W, length(j), lower.tail = FALSE)) }
w("Within-group contrast, beta1 (SE): ")
for (g in names(ct)) { j <- match(paste0("beta1:", Y), nm)
  w("- ", g, ": ", paste(sprintf("%+.3f (%.3f)", ct[[g]]$d[j], sqrt(diag(ct[[g]]$V))[j]), collapse = ", "),
    " | block p ", fp(wald(ct[[g]]$d, ct[[g]]$V, j)["p"])) }
w(""); w("| comparison | block | W (4 df) | p |"); w("|---|---|---|---|")
for (g in c("Recent/mid-term", "Established", "Immigrant (all)")) {
  D <- ct[[g]]$d - ct[["Non-immigrant"]]$d; VD <- ct[[g]]$V + ct[["Non-immigrant"]]$V
  for (k in CMP) { t <- wald(D, VD, match(paste0(k, ":", Y), nm))
    w("| ", g, " vs Non-immigrant | ", k, " | ", sprintf("%.2f", t["W"]), " | ", fp(t["p"]), " |") }
  j <- match(paste0("beta1:", Y), nm)
  w("|  | beta1 elementwise | ", paste(sprintf("%+.3f (%.3f)", D[j], sqrt(diag(VD))[j]), collapse = "; "), " | ",
    paste(fp(2 * pnorm(-abs(D[j] / sqrt(diag(VD))[j]))), collapse = "/"), " |")
}
w(""); w("Stratum sizes are not written here (secondary-disclosure control).")

# ===========================================================================
w(""); w("## PART 4: instability tests at inner nodes of the full-sample slope tree"); w("")
mob <- readRDS(file.path(OUT, "_mob_final_model.rds"))
tr  <- mob[["D1 slope"]]$tree
for (i in setdiff(nodeids(tr), nodeids(tr, terminal = TRUE))) {
  nd <- nodeapply(tr, i, identity)[[1]]; ti <- info_node(nd)
  sp <- split_node(nd); vn <- names(tr$data)[sp$varid]
  o <- order(ti$test["p.value", ])
  w("- node ", i, " (split on ", vn, ", ", sup(ti$nobs %/% 3), " persons): ",
    paste(sprintf("%s p=%.3g", colnames(ti$test)[o], ti$test["p.value", o]), collapse = ", "))
}
close(con)
cat("\nwritten: DataPrep/out/MANUSCRIPT_REVISIONS.md\n")
