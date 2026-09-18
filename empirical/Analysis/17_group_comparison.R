# ---------------------------------------------------------------------------
# 17_group_comparison.R
#
# The final foundational ERA-GEE model fitted SEPARATELY in the two groups
# defined by the CLSA Immigration Flag (SDC_FIMM_COM), with a formal
# comparison of every parameter between groups.
#
# MODEL (identical to Analysis/15, fitted within each group)
#     E(Y_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c)
#
#   t    = data-collection wave index (0, 1, 2)
#   A_c  = age at the CLSA baseline visit, centred at the FIXED reference age
#          of the balanced sample (62.15 years) in BOTH groups, so that
#          beta_1q and beta_2q refer to the same age in each group. Centring
#          each group at its own mean would compare quantities defined at
#          different ages.
#
# Fitting the model separately in the two groups is equivalent to a single
# model in which every component is fully interacted with immigrant status.
#
# COMPARISON. The two groups are disjoint sets of participants and clustering
# is within participant, so the two estimators are independent:
#     Var(theta_imm - theta_non) = Var(theta_imm) + Var(theta_non).
# Differences are tested elementwise and with joint Wald tests: one per
# parameter block (4 df, across the four outcomes) and one overall (16 df).
#
# Run from the project root:
#     Rscript Analysis/17_group_comparison.R
# Writes (local only): DataPrep/out/GROUP_COMPARISON.md
# ---------------------------------------------------------------------------

source("Analysis/era_gee_components.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c(y_walk = "Walking", y_lsport = "Light", y_msport = "Moderate", y_ssport = "Strenuous")
CMP <- c("beta1", "beta2", "gamma", "delta")

bal  <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
alls <- readRDS(file.path(OUT, "analytic_long_all.rds"))
raw  <- readRDS(file.path(OUT, "immigration_items_raw.rds"))

AGE_C <- mean(bal$z_age[bal$wave == 0])          # fixed reference age, both groups

add_flag <- function(d) {
  f <- raw$SDC_FIMM_COM[match(d$entity_id, raw$entity_id)]
  d$grp <- factor(ifelse(f == "1", "Immigrant", ifelse(f == "2", "Non-immigrant", NA)),
                  levels = c("Non-immigrant", "Immigrant"))
  d[!is.na(d$grp), ]
}
bal <- add_flag(bal); alls <- add_flag(alls)

design <- function(d) {
  Ac <- d$z_age - AGE_C
  cbind(beta1 = 1, beta2 = d$X_time, gamma = Ac, delta = d$X_time * Ac)
}

fit_groups <- function(d) {
  lapply(split(seq_len(nrow(d)), d$grp), function(ix) {
    dd <- d[ix, ]
    era_gee_components(dd[, Y], design(dd), dd$id)
  })
}

# --- difference of two independent parameter vectors ------------------------
compare <- function(f_non, f_imm) {
  nm <- rownames(f_non$V)
  dif <- f_imm$theta - f_non$theta
  Vd  <- f_imm$V + f_non$V
  se  <- sqrt(diag(Vd))
  z   <- dif / se
  wald <- function(which) {
    i <- match(which, nm)
    W <- as.numeric(t(dif[i]) %*% solve(Vd[i, i, drop = FALSE]) %*% dif[i])
    c(W = W, df = length(i), p = pchisq(W, length(i), lower.tail = FALSE))
  }
  list(nm = nm, dif = dif, se = se, z = z, p = 2 * pnorm(-abs(z)),
       block = sapply(CMP, function(k) wald(paste0(k, ":", Y))),
       overall = wald(nm))
}

# ---------------------------------------------------------------------------
fits <- fit_groups(bal)
cmp  <- compare(fits[["Non-immigrant"]], fits[["Immigrant"]])
fitsA <- fit_groups(alls)
cmpA  <- compare(fitsA[["Non-immigrant"]], fitsA[["Immigrant"]])

# --- model-implied per-wave change by age, within group --------------------
per_wave <- function(f, off) sapply(Y, function(q) f$B["beta2", q] + off * f$B["delta", q])
lvl_wave <- function(f, off) sapply(Y, function(q) f$B["beta1", q] + off * f$B["gamma", q])

# ---------------------------------------------------------------------------
fp  <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
b0  <- bal[bal$wave == 0, ]

con <- file(file.path(OUT, "GROUP_COMPARISON.md"), open = "wt")
w <- function(...) writeLines(paste0(...), con)

w("# Foundational ERA-GEE model fitted separately by immigrant status")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `Analysis/17_group_comparison.R`")
w("")
w("Groups from the CLSA Immigration Flag (`SDC_FIMM_COM`). Both groups use the")
w("same fixed reference age of ", sprintf("%.2f", AGE_C), " years, so beta1 and beta2 refer to the")
w("same age in each group. Estimates are in scale points on the 1-4 response scale.")
w("")
w("## Groups (balanced sample)")
w("")
w("| | Non-immigrant | Immigrant |"); w("|---|---|---|")
g <- split(b0, b0$grp)
w("| Participants | ", format(nrow(g[["Non-immigrant"]]), big.mark = ","), " | ",
  format(nrow(g[["Immigrant"]]), big.mark = ","), " |")
w("| Person-waves | ", format(3 * nrow(g[["Non-immigrant"]]), big.mark = ","), " | ",
  format(3 * nrow(g[["Immigrant"]]), big.mark = ","), " |")
w("| Baseline-visit age, M (SD) | ", sprintf("%.2f (%.2f)", mean(g[[1]]$z_age), sd(g[[1]]$z_age)),
  " | ", sprintf("%.2f (%.2f)", mean(g[[2]]$z_age), sd(g[[2]]$z_age)), " |")
w("| Female, % | ", sprintf("%.1f", 100 * mean(g[[1]]$z_sex == "Female")), " | ",
  sprintf("%.1f", 100 * mean(g[[2]]$z_sex == "Female")), " |")
w("| Racialized, % | ", sprintf("%.1f", 100 * mean(g[[1]]$z_racialized == "Racialized", na.rm = TRUE)),
  " | ", sprintf("%.1f", 100 * mean(g[[2]]$z_racialized == "Racialized", na.rm = TRUE)), " |")
w("")
w("## Observed wave means by group")
w("")
w("| Outcome | ", paste(paste0("Non W", 1:3), collapse = " | "), " | ",
  paste(paste0("Imm W", 1:3), collapse = " | "), " |")
w("|---|---|---|---|---|---|---|")
for (q in Y) {
  m <- sapply(c("Non-immigrant", "Immigrant"), function(gg)
    sapply(0:2, function(ww) mean(bal[[q]][bal$grp == gg & bal$wave == ww])))
  w("| ", YL[[q]], " | ", paste(sprintf("%.3f", m[, 1]), collapse = " | "), " | ",
    paste(sprintf("%.3f", m[, 2]), collapse = " | "), " |")
}
w("")
w("## Parameter estimates and between-group differences")
w("")
w("| Parameter | Outcome | Non-immigrant | Immigrant | Difference (Imm - Non) | 95% CI | p |")
w("|---|---|---|---|---|---|---|")
for (k in CMP) for (q in Y) {
  i  <- match(paste0(k, ":", q), cmp$nm)
  dg <- if (k %in% c("beta1", "beta2")) 4 else 5
  f1 <- sprintf(paste0("%+.", dg, "f (%.", dg, "f)"), fits[[1]]$B[k, q], fits[[1]]$SE[k, q])
  f2 <- sprintf(paste0("%+.", dg, "f (%.", dg, "f)"), fits[[2]]$B[k, q], fits[[2]]$SE[k, q])
  w("| ", k, " | ", YL[[q]], " | ", f1, " | ", f2, " | ",
    sprintf(paste0("%+.", dg, "f (%.", dg, "f)"), cmp$dif[i], cmp$se[i]), " | [",
    sprintf(paste0("%+.", dg, "f"), cmp$dif[i] - 1.96 * cmp$se[i]), ", ",
    sprintf(paste0("%+.", dg, "f"), cmp$dif[i] + 1.96 * cmp$se[i]), "] | ", fp(cmp$p[i]), " |")
}
w("")
w("Standard errors in parentheses. Individual p values are unadjusted for the")
w("16 comparisons; the joint tests below are the primary evidence.")
w("")
w("## Joint tests of between-group equality")
w("")
w("| Parameter block | W | df | p |"); w("|---|---|---|---|")
for (k in CMP) w("| ", k, " | ", sprintf("%.2f", cmp$block["W", k]), " | 4 | ", fp(cmp$block["p", k]), " |")
w("| all 16 parameters | ", sprintf("%.2f", cmp$overall["W"]), " | 16 | ", fp(cmp$overall["p"]), " |")
w("")
w("## Model-implied per-wave change by baseline-visit age and group")
w("")
w("| Baseline-visit age | Group | ", paste(YL[Y], collapse = " | "), " |")
w("|---|---|---|---|---|---|")
for (off in c(-10, 0, 10)) for (gi in 1:2)
  w("| ", sprintf("%.2f", AGE_C + off), " | ", names(fits)[gi], " | ",
    paste(sprintf("%+.4f", per_wave(fits[[gi]], off)), collapse = " | "), " |")
w("")
w("## Sensitivity: all available complete person-wave records")
w("")
w("| Parameter block | W | df | p |"); w("|---|---|---|---|")
for (k in CMP) w("| ", k, " | ", sprintf("%.2f", cmpA$block["W", k]), " | 4 | ", fp(cmpA$block["p", k]), " |")
w("| all 16 parameters | ", sprintf("%.2f", cmpA$overall["W"]), " | 16 | ", fp(cmpA$overall["p"]), " |")
w("")
close(con)

# --- console ---------------------------------------------------------------
cat("reference age:", round(AGE_C, 2), "\n")
cat("\ngroup sizes (balanced):\n"); print(table(b0$grp))
cat("\nNon-immigrant B:\n"); print(round(fits[["Non-immigrant"]]$B, 5))
cat("Non-immigrant SE:\n"); print(round(fits[["Non-immigrant"]]$SE, 5))
cat("\nImmigrant B:\n"); print(round(fits[["Immigrant"]]$B, 5))
cat("Immigrant SE:\n"); print(round(fits[["Immigrant"]]$SE, 5))
cat("\ndifferences (Imm - Non) with SE and p:\n")
print(data.frame(parameter = cmp$nm, diff = round(cmp$dif, 5),
                 se = round(cmp$se, 5), p = signif(cmp$p, 3)), row.names = FALSE)
cat("\njoint tests (balanced):\n"); print(round(cmp$block, 4)); print(round(cmp$overall, 4))
cat("\njoint tests (all available):\n"); print(round(cmpA$block, 4)); print(round(cmpA$overall, 4))
cat("\nper-wave change by age and group:\n")
for (off in c(-10, 0, 10)) for (gi in 1:2)
  cat(sprintf("  age %.2f  %-14s %s\n", AGE_C + off, names(fits)[gi],
              paste(sprintf("%+.4f", per_wave(fits[[gi]], off)), collapse = "  ")))
cat("\nfitted Wave-1 level by age and group:\n")
for (off in c(-10, 0, 10)) for (gi in 1:2)
  cat(sprintf("  age %.2f  %-14s %s\n", AGE_C + off, names(fits)[gi],
              paste(sprintf("%.4f", lvl_wave(fits[[gi]], off)), collapse = "  ")))
cat("\nwritten: DataPrep/out/GROUP_COMPARISON.md\n")
