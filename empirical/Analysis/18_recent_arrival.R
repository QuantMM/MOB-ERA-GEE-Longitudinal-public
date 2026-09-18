# ---------------------------------------------------------------------------
# 18_recent_arrival.R
#
# The final foundational ERA-GEE model fitted separately by SETTLEMENT
# DURATION, with the focus on participants who had been in Canada fewer than
# 10 years at the CLSA baseline visit.
#
# MODEL (identical to Analysis/15 and 17, fitted within each group)
#     E(Y_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c)
#     A_c centred at the FIXED reference age of 62.15 years in every group.
#
# GROUPS
#     Non-immigrant                (CLSA Immigration Flag = 2)
#     Immigrant, < 10 years        (SDC_DRES_COM < 10)      <- focus
#     Immigrant, 10-19 years
#     Immigrant, 20-39 years
#     Immigrant, 40+ years
#
# WHY THIS IS A POWER QUESTION FIRST. The < 10 year group is very small, so
# the comparison is reported together with the minimum difference the design
# could detect: for two independent groups the 80%-power two-sided detectable
# difference at alpha = .05 is approximately 2.8 x SE(difference). A null
# result here means "not detectable", not "not present".
#
# EXTRAPOLATION. Recent arrivals are much younger than the cohort mean, so
# beta_1 and beta_2 evaluated at 62.15 years may sit in the upper tail of this
# group's age range. The age distribution of each group is reported, and
# model-implied quantities are also given at 52.15 years, where every group
# has support.
#
# Run from the project root:
#     Rscript Analysis/18_recent_arrival.R
# Writes (local only): DataPrep/out/RECENT_ARRIVAL.md
# ---------------------------------------------------------------------------

source("Analysis/era_gee_components.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c(y_walk = "Walking", y_lsport = "Light", y_msport = "Moderate", y_ssport = "Strenuous")
CMP <- c("beta1", "beta2", "gamma", "delta")

bal  <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
alls <- readRDS(file.path(OUT, "analytic_long_all.rds"))
raw  <- readRDS(file.path(OUT, "immigration_items_raw.rds"))
AGE_C <- mean(bal$z_age[bal$wave == 0])

LV <- c("Non-immigrant", "< 10 yr", "10-19 yr", "20-39 yr", "40+ yr")
add_grp <- function(d) {
  i  <- match(d$entity_id, raw$entity_id)
  fl <- raw$SDC_FIMM_COM[i]
  yr <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[i]))
  g  <- ifelse(fl == "2", "Non-immigrant",
        ifelse(fl == "1" & !is.na(yr) & yr < 996,
               ifelse(yr < 10, "< 10 yr",
               ifelse(yr < 20, "10-19 yr",
               ifelse(yr < 40, "20-39 yr", "40+ yr"))), NA))
  d$grp <- factor(g, levels = LV)
  d[!is.na(d$grp), ]
}
bal <- add_grp(bal); alls <- add_grp(alls)

design <- function(d) {
  Ac <- d$z_age - AGE_C
  cbind(beta1 = 1, beta2 = d$X_time, gamma = Ac, delta = d$X_time * Ac)
}
fit_one <- function(d) era_gee_components(d[, Y], design(d), d$id)

fit_groups <- function(d, groups) {
  out <- list()
  for (g in groups) { dd <- d[d$grp %in% g, ]; out[[paste(g, collapse = "+")]] <- fit_one(dd) }
  out
}

# --- difference between two independent groups -----------------------------
compare <- function(fa, fb) {          # fb - fa
  nm  <- rownames(fa$V)
  dif <- fb$theta - fa$theta
  Vd  <- fa$V + fb$V
  se  <- sqrt(diag(Vd))
  wald <- function(which) { i <- match(which, nm)
    W <- as.numeric(t(dif[i]) %*% solve(Vd[i, i, drop = FALSE]) %*% dif[i])
    c(W = W, df = length(i), p = pchisq(W, length(i), lower.tail = FALSE)) }
  list(nm = nm, dif = dif, se = se, p = 2 * pnorm(-abs(dif / se)),
       mdd = 2.8 * se,                                  # 80% power, alpha .05
       block = sapply(CMP, function(k) wald(paste0(k, ":", Y))),
       overall = wald(nm))
}

fits  <- fit_groups(bal, as.list(LV))
fits[["Immigrant 10+ yr"]] <- fit_one(bal[bal$grp %in% LV[3:5], ])
cmp_non  <- compare(fits[["Non-immigrant"]],    fits[["< 10 yr"]])
cmp_long <- compare(fits[["Immigrant 10+ yr"]], fits[["< 10 yr"]])

fitsA <- fit_groups(alls, as.list(LV))
cmpA  <- compare(fitsA[["Non-immigrant"]], fitsA[["< 10 yr"]])

fp <- function(p) ifelse(p < .001, "<.001", sub("^0", "", sprintf("%.3f", p)))
b0 <- bal[bal$wave == 0, ]

con <- file(file.path(OUT, "RECENT_ARRIVAL.md"), open = "wt")
w <- function(...) writeLines(paste0(...), con)

w("# Final ERA-GEE model by settlement duration, focusing on arrivals within 10 years")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `Analysis/18_recent_arrival.R`")
w("")
w("All groups use the same fixed reference age of ", sprintf("%.2f", AGE_C), " years. Estimates are in")
w("scale points on the 1-4 response scale.")
w("")
w("## Groups (balanced sample)")
w("")
w("| Group | n | Age M (SD) | Age range | % above 62.15 | % female | % racialized |")
w("|---|---|---|---|---|---|---|")
for (g in LV) { s <- b0[b0$grp == g, ]
  w("| ", g, " | ", format(nrow(s), big.mark = ","), " | ",
    sprintf("%.1f (%.1f)", mean(s$z_age), sd(s$z_age)), " | ",
    sprintf("%d-%d", min(s$z_age), max(s$z_age)), " | ",
    sprintf("%.0f", 100 * mean(s$z_age > AGE_C)), " | ",
    sprintf("%.0f", 100 * mean(s$z_sex == "Female")), " | ",
    sprintf("%.0f", 100 * mean(s$z_racialized == "Racialized", na.rm = TRUE)), " |") }
w("")
w("## Observed wave means")
w("")
w("| Group | ", paste(rep(YL[Y], each = 1), collapse = " | "), " |")
w("|---|---|---|---|---|")
for (g in LV) w("| ", g, " | ", paste(sapply(Y, function(q)
  paste(sprintf("%.2f", sapply(0:2, function(ww) mean(bal[[q]][bal$grp == g & bal$wave == ww]))),
        collapse = "/")), collapse = " | "), " |")
w("")
w("Cells give the Wave 1 / Wave 2 / Wave 3 means.")
w("")
w("## Estimates by group")
w("")
for (k in CMP) {
  dg <- if (k %in% c("beta1", "beta2")) 3 else 5
  w("**", k, "**")
  w("")
  w("| Group | ", paste(YL[Y], collapse = " | "), " |"); w("|---|---|---|---|---|")
  for (g in c(LV, "Immigrant 10+ yr")) {
    f <- fits[[g]]
    w("| ", g, " | ", paste(sprintf(paste0("%+.", dg, "f (%.", dg, "f)"),
      f$B[k, Y], f$SE[k, Y]), collapse = " | "), " |")
  }
  w("")
}
w("## Comparison: arrivals within 10 years vs each reference group")
w("")
for (nmc in c("vs non-immigrants", "vs immigrants settled 10+ years")) {
  cc <- if (nmc == "vs non-immigrants") cmp_non else cmp_long
  w("### ", nmc)
  w("")
  w("| Parameter block | W | df | p |"); w("|---|---|---|---|")
  for (k in CMP) w("| ", k, " | ", sprintf("%.2f", cc$block["W", k]), " | 4 | ", fp(cc$block["p", k]), " |")
  w("| all 16 | ", sprintf("%.2f", cc$overall["W"]), " | 16 | ", fp(cc$overall["p"]), " |")
  w("")
  w("| Parameter | Outcome | Difference (SE) | p | Detectable difference |")
  w("|---|---|---|---|---|")
  for (k in c("beta1", "beta2")) for (q in Y) { i <- match(paste0(k, ":", q), cc$nm)
    w("| ", k, " | ", YL[[q]], " | ", sprintf("%+.3f (%.3f)", cc$dif[i], cc$se[i]), " | ",
      fp(cc$p[i]), " | ", sprintf("%.3f", cc$mdd[i]), " |") }
  w("")
}
w("The last column is the difference this design could detect with 80% power at")
w("alpha = .05, approximately 2.8 times the standard error of the difference.")
w("")
w("## Model-implied per-wave change at two reference ages")
w("")
for (off in c(-10, 0)) {
  w("**Baseline-visit age ", sprintf("%.2f", AGE_C + off), "**")
  w("")
  w("| Group | ", paste(YL[Y], collapse = " | "), " |"); w("|---|---|---|---|---|")
  for (g in LV) { f <- fits[[g]]
    w("| ", g, " | ", paste(sprintf("%+.3f", f$B["beta2", Y] + off * f$B["delta", Y]),
                            collapse = " | "), " |") }
  w("")
}
w("## Sensitivity: all available complete person-wave records")
w("")
w("Arrivals within 10 years: n = ", length(unique(alls$id[alls$grp == "< 10 yr"])), ".")
w("")
w("| Parameter block | W | df | p |"); w("|---|---|---|---|")
for (k in CMP) w("| ", k, " | ", sprintf("%.2f", cmpA$block["W", k]), " | 4 | ", fp(cmpA$block["p", k]), " |")
w("| all 16 | ", sprintf("%.2f", cmpA$overall["W"]), " | 16 | ", fp(cmpA$overall["p"]), " |")
w("")
close(con)

# --- console ---------------------------------------------------------------
cat("group sizes (balanced):\n"); print(table(b0$grp))
cat("\nage by group:\n")
for (g in LV) { s <- b0[b0$grp == g, ]
  cat(sprintf("  %-14s n=%6d  age %.1f (%.1f)  range %d-%d  %.0f%% above %.2f  %.0f%%F  %.0f%%rac\n",
      g, nrow(s), mean(s$z_age), sd(s$z_age), min(s$z_age), max(s$z_age),
      100 * mean(s$z_age > AGE_C), AGE_C, 100 * mean(s$z_sex == "Female"),
      100 * mean(s$z_racialized == "Racialized", na.rm = TRUE))) }
cat("\nobserved wave means, < 10 yr group:\n")
for (q in Y) cat(sprintf("  %-9s %s\n", YL[[q]],
  paste(sprintf("%.3f", sapply(0:2, function(ww) mean(bal[[q]][bal$grp == "< 10 yr" & bal$wave == ww]))), collapse = "  ")))
for (g in c(LV, "Immigrant 10+ yr")) {
  cat("\n===== ", g, " =====\n"); print(round(fits[[g]]$B, 5)); cat("SE:\n"); print(round(fits[[g]]$SE, 5)) }
cat("\n--- < 10 yr vs non-immigrants ---\n"); print(round(cmp_non$block, 4)); print(round(cmp_non$overall, 4))
print(data.frame(p = cmp_non$nm, diff = round(cmp_non$dif, 4), se = round(cmp_non$se, 4),
                 pval = signif(cmp_non$p, 3), mdd = round(cmp_non$mdd, 3)), row.names = FALSE)
cat("\n--- < 10 yr vs immigrants 10+ yr ---\n"); print(round(cmp_long$block, 4)); print(round(cmp_long$overall, 4))
cat("\n--- sensitivity (all available) ---\n"); print(round(cmpA$block, 4)); print(round(cmpA$overall, 4))
cat("\nper-wave change at 52.15 and 62.15:\n")
for (off in c(-10, 0)) for (g in LV) cat(sprintf("  age %.2f  %-14s %s\n", AGE_C + off, g,
  paste(sprintf("%+.4f", fits[[g]]$B["beta2", Y] + off * fits[[g]]$B["delta", Y]), collapse = "  ")))
cat("\nwritten: DataPrep/out/RECENT_ARRIVAL.md\n")
