# ---------------------------------------------------------------------------
# 19_three_group.R
#
# FINAL ANALYSIS. Three parts, all built on the final foundational ERA-GEE
# model  E(Y_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c),
# with A_c centred at the fixed reference age of 62.15 years everywhere.
#
#   GROUPS   Non-immigrant            (CLSA Immigration Flag = 2)
#            Recent/mid-term          (Flag = 1, SDC_DRES_COM  < 20 years)
#            Established/long-term    (Flag = 1, SDC_DRES_COM >= 20 years)
#
#   PART A   The model fitted separately in each group; all pairwise parameter
#            comparisons and a three-group omnibus test.
#   PART B   MOB within each group, with an IDENTICAL candidate set Z so the
#            groups are comparable: sex, racialization, education, income
#            adequacy, dwelling ownership. gamma and delta stay global within
#            each group; beta_1 and beta_2 are node-specific and tested.
#   PART C   Prespecified stratification. The same two contrasts (sex,
#            racialization) are estimated in every group whether or not MOB
#            selects them, so the cross-group comparison does not depend on
#            tree shape, which is confounded with group size.
#
# WHY PART C EXISTS. Recent/mid-term immigrants number a few hundred and
# non-immigrants twenty thousand. With minsize = 100 persons the smaller group
# can support only a couple of splits, so "fewer splits" is partly a statement
# about sample size. Root instability and prespecified contrasts are the
# power-transparent comparisons; tree shape is not.
#
# Trees are grown permissively and ARE NOT PRUNED. Counts below 6 suppressed.
#
# Run from the project root:
#     Rscript Analysis/19_three_group.R
# Writes (local only): DataPrep/out/THREE_GROUP.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("FINAL_MOB-ERA-GEE_Code/Empirical Application/FINAL_HELPERS_MOB_ERA_GEE.R")
source("Analysis/ERA_GEE_fixed.R")
source("Analysis/nuisance_adjusted_slope.R")
source("Analysis/palm_era_gee.R")
source("Analysis/era_gee_components.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c(y_walk = "Walking", y_lsport = "Light", y_msport = "Moderate", y_ssport = "Strenuous")
CMP <- c("beta1", "beta2", "gamma", "delta")
Q   <- length(Y); FAM <- rep("gaussian", Q)
Z   <- c("z_sex", "z_racialized", "z_educ4", "z_incneeds", "z_own")
MINSIZE <- 100L; MAXDEPTH <- 4L; ALPHA <- 0.05
GRP <- c("Non-immigrant", "Recent/mid-term", "Established")

bal <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
raw <- readRDS(file.path(OUT, "immigration_items_raw.rds"))
AGE_C <- mean(bal$z_age[bal$wave == 0])

i  <- match(bal$entity_id, raw$entity_id)
fl <- raw$SDC_FIMM_COM[i]; yr <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[i]))
bal$grp <- factor(ifelse(fl == "2", "Non-immigrant",
                  ifelse(fl == "1" & !is.na(yr) & yr < 996,
                         ifelse(yr < 20, "Recent/mid-term", "Established"), NA)), levels = GRP)
bal <- bal[!is.na(bal$grp), ]

sup  <- function(n) if (n > 0 && n < 6) "<6" else format(n, big.mark = ",", trim = TRUE)
fp   <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
line <- function(ch = "-") cat(strrep(ch, 78), "\n")
design <- function(d) { Ac <- d$z_age - AGE_C
  cbind(beta1 = 1, beta2 = d$X_time, gamma = Ac, delta = d$X_time * Ac) }
fit_one <- function(d) era_gee_components(d[, Y], design(d), d$id)

# ===========================================================================
cat("\n"); line("#"); cat("PART A -- three-group parameter comparison\n"); line("#")
fits <- lapply(setNames(GRP, GRP), function(g) fit_one(bal[bal$grp == g, ]))
b0 <- bal[bal$wave == 0, ]
for (g in GRP) { s <- b0[b0$grp == g, ]
  cat(sprintf("  %-16s n=%6s  age %.1f (%.1f)  %.0f%%F  %.0f%%rac\n", g, sup(nrow(s)),
      mean(s$z_age), sd(s$z_age), 100*mean(s$z_sex == "Female"),
      100*mean(s$z_racialized == "Racialized", na.rm = TRUE))) }

cat("\nobserved wave means (W1/W2/W3):\n")
for (g in GRP) { cat(sprintf("  %-16s %s\n", g, paste(sapply(Y, function(q)
  paste(sprintf("%.3f", sapply(0:2, function(w) mean(bal[[q]][bal$grp == g & bal$wave == w]))),
        collapse = "/")), collapse = "  "))) }

for (g in GRP) { cat("\n=====", g, "=====\n"); print(round(fits[[g]]$B, 5))
  cat("SE:\n"); print(round(fits[[g]]$SE, 5)) }

pair <- function(fa, fb) {                       # fb - fa
  nm <- rownames(fa$V); dif <- fb$theta - fa$theta; Vd <- fa$V + fb$V
  se <- sqrt(diag(Vd))
  wald <- function(wh) { j <- match(wh, nm)
    W <- as.numeric(t(dif[j]) %*% solve(Vd[j, j, drop = FALSE]) %*% dif[j])
    c(W = W, df = length(j), p = pchisq(W, length(j), lower.tail = FALSE)) }
  list(nm = nm, dif = dif, se = se, p = 2*pnorm(-abs(dif/se)), mdd = 2.8*se,
       block = sapply(CMP, function(k) wald(paste0(k, ":", Y))), overall = wald(nm))
}
PAIRS <- list(`Recent vs Non` = c(1, 2), `Established vs Non` = c(1, 3),
              `Recent vs Established` = c(3, 2))
cmp <- lapply(PAIRS, function(k) pair(fits[[GRP[k[1]]]], fits[[GRP[k[2]]]]))
for (nmp in names(cmp)) { cat("\n---", nmp, "---\n")
  print(round(cmp[[nmp]]$block, 4)); print(round(cmp[[nmp]]$overall, 4)) }

# three-group omnibus: both contrasts against non-immigrants, jointly
d12 <- fits[[2]]$theta - fits[[1]]$theta; d13 <- fits[[3]]$theta - fits[[1]]$theta
V1 <- fits[[1]]$V
Vom <- rbind(cbind(fits[[2]]$V + V1, V1), cbind(V1, fits[[3]]$V + V1))
dom <- c(d12, d13); k <- length(fits[[1]]$theta)
omni <- function(idx) { j <- c(idx, idx + k)
  W <- as.numeric(t(dom[j]) %*% solve(Vom[j, j, drop = FALSE]) %*% dom[j])
  c(W = W, df = length(j), p = pchisq(W, length(j), lower.tail = FALSE)) }
cat("\n--- three-group omnibus (all parameters equal across the three groups) ---\n")
ob <- sapply(CMP, function(cc) omni(match(paste0(cc, ":", Y), rownames(V1))))
print(round(ob, 4)); print(round(omni(seq_len(k)), 4))

cat("\nminimum detectable difference, Recent vs Non (beta2):\n")
mdd <- cmp[["Recent vs Non"]]$mdd[match(paste0("beta2:", Y), cmp[[1]]$nm)]
cat(sprintf("  %-10s %.3f\n", YL[Y], mdd), sep = "")

# ===========================================================================
cat("\n"); line("#"); cat("PART B -- MOB within each group (identical Z)\n"); line("#")
ctrl <- function(prm) mob_control(verbose = FALSE, maxdepth = MAXDEPTH,
                                  minsize = MINSIZE * 3L, alpha = ALPHA, parm = prm)
split_desc <- function(tr) {
  inner <- setdiff(nodeids(tr), nodeids(tr, terminal = TRUE))
  if (!length(inner)) return("(no split)")
  paste(sapply(inner, function(i) { nd <- nodeapply(tr, i, identity)[[1]]; sp <- split_node(nd)
    vn <- names(tr$data)[sp$varid]; zz <- tr$data[[sp$varid]]
    kk <- sp$index
    sprintf("node %d: %s {%s}|{%s}", i, vn,
            paste(levels(zz)[which(kk == 1L)], collapse = ","),
            paste(levels(zz)[which(kk == 2L)], collapse = ","))}), collapse = " || ")
}
mobres <- list()
for (g in GRP) {
  dd <- droplevels(bal[bal$grp == g, ]); dd$rowid <- seq_len(nrow(dd))
  np <- length(unique(dd$id))
  cat("\n"); line("="); cat(sprintf("%s : %s persons  (minsize = %d persons)\n", g, sup(np), MINSIZE)); line("=")
  for (dg in c("joint", "slope", "intercept")) {
    prm <- switch(dg, joint = 1:(2*Q), intercept = 1:Q, slope = NULL)
    f  <- palm_era_gee_mob(dd, Z, dg, Y, "X_time", FAM, age_var = "z_age",
                           age_center = AGE_C, control = ctrl(prm))
    ti <- nodeapply(f$tree, 1, function(n) info_node(n))[[1]]
    o  <- order(ti$test["p.value", ])
    cat(sprintf("\n  %s | %s: %d terminal nodes | %s\n", g, dg, width(f$tree), split_desc(f$tree)))
    cat("    root instability: ")
    cat(paste(sprintf("%s=%.3g", colnames(ti$test)[o], ti$test["p.value", o]), collapse = "  "), "\n")
    if (width(f$tree) > 1L) {
      node <- f$node; off <- f$offset
      for (nd in sort(unique(node[!is.na(node)]))) {
        ix <- which(!is.na(node) & node == nd); s <- dd[ix, ]
        n2 <- length(unique(s$id)); if (n2 < 6) { cat("     node", nd, "[suppressed]\n"); next }
        m <- ERA_GEE_fixed(as.matrix(s[, Y]) - off[ix, , drop = FALSE], s$X_time, FAM,
                           "independence", y_scale = "raw")
        cat(sprintf("     node %-2d n=%6s %3.0f%%F %3.0f%%rac  b2 %s\n", nd, sup(n2),
                    100*mean(s$z_sex == "Female"), 100*mean(s$z_racialized == "Racialized", na.rm = TRUE),
                    paste(sprintf("%+7.4f", m$B["slope", ]), collapse = " ")))
      }
    }
    mobres[[paste(g, dg)]] <- f
  }
}

# ===========================================================================
cat("\n"); line("#"); cat("PART C -- prespecified stratification (same contrasts in every group)\n"); line("#")
strat <- list(sex = c("z_sex", "Female", "Male"),
              racialization = c("z_racialized", "White", "Racialized"))
ctab <- list()
for (sname in names(strat)) {
  v <- strat[[sname]][1]; l1 <- strat[[sname]][2]; l2 <- strat[[sname]][3]
  cat("\n=== contrast:", l2, "minus", l1, "(", sname, ") ===\n")
  for (g in GRP) {
    dg1 <- bal[bal$grp == g & !is.na(bal[[v]]) & bal[[v]] == l1, ]
    dg2 <- bal[bal$grp == g & !is.na(bal[[v]]) & bal[[v]] == l2, ]
    n1 <- length(unique(dg1$id)); n2 <- length(unique(dg2$id))
    if (min(n1, n2) < 30) { cat(sprintf("  %-16s skipped (n = %s / %s)\n", g, sup(n1), sup(n2))); next }
    f1 <- fit_one(dg1); f2 <- fit_one(dg2); cc <- pair(f1, f2)
    ctab[[paste(sname, g)]] <- list(cmp = cc, n1 = n1, n2 = n2)
    cat(sprintf("  %-16s n = %s vs %s | block p: %s\n", g, sup(n1), sup(n2),
                paste(sprintf("%s=%s", CMP, sapply(CMP, function(k) fp(cc$block["p", k]))), collapse = "  ")))
    j <- match(paste0("beta2:", Y), cc$nm)
    cat(sprintf("     beta2 contrast %s\n",
                paste(sprintf("%+.4f(%.4f)", cc$dif[j], cc$se[j]), collapse = " ")))
  }
}

saveRDS(list(fits = fits, cmp = cmp, omnibus = list(block = ob, all = omni(seq_len(k))),
             mob = mobres, strat = ctab), file.path(OUT, "_three_group.rds"))
cat("\nsaved: DataPrep/out/_three_group.rds\n")
