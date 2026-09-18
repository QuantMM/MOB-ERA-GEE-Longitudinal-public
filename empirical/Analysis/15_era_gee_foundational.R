# ---------------------------------------------------------------------------
# 15_era_gee_foundational.R
#
# The foundational ERA-GEE model on the CLSA data, WITHOUT any partitioning.
#
# DEFINITIONS (agreed wording)
#   t    = 0, 1, 2 is the DATA-COLLECTION WAVE INDEX, not elapsed calendar
#          time. In the CLSA the two between-wave intervals differ (median
#          about 1.5 and 2.8 years), so beta_2q is the average linear change
#          per wave, not an annual change, and a single linear summary of the
#          pattern across the three waves rather than a constant true slope.
#   A_c  = age at the CLSA baseline visit, centred at the analytic-sample
#          mean (person level). This is NOT age at the first physical-activity
#          measurement: Wave 1 PA comes from a later telephone interview.
#
# MODELS (Gaussian, identity link, working independence, all component
# weights fixed at 1, person-clustered sandwich SEs; Q = 4 outcomes)
#   M2  g(mu_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c)
#       MAIN age-adjusted model: evaluates age-related heterogeneity in
#       longitudinal change (delta_q).
#   M0  g(mu_itq) = beta_1q + beta_2q t               (reference)
#   M1  g(mu_itq) = beta_1q + beta_2q t + gamma_q A_c (reference; the
#       comparison with M2 documents why the age x wave term is needed)
#
# Samples: balanced (complete at all three waves; primary) and all available
# person-waves (sensitivity).
#
# Run from the project root:
#     Rscript Analysis/15_era_gee_foundational.R
# Writes (local only): DataPrep/out/ERA_GEE_FOUNDATIONAL.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(sandwich))
source("FINAL_MOB-ERA-GEE_Code/Empirical Application/FINAL_HELPERS_MOB_ERA_GEE.R")
source("Analysis/ERA_GEE_fixed.R")
source("Analysis/era_gee_components.R")

OUT  <- "DataPrep/out"
DATA <- "CLSA Datafiles for analysis"
Y    <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL   <- c(y_walk = "Walking", y_lsport = "Light", y_msport = "Moderate", y_ssport = "Strenuous")

bal  <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
alls <- readRDS(file.path(OUT, "analytic_long_all.rds"))

# ONE FIXED REFERENCE AGE FOR BOTH SAMPLES.
# With the age x wave component in the model, beta_2q is the per-wave change AT
# THE REFERENCE AGE. Centring each sample at its own mean would therefore make
# the sensitivity comparison a comparison of slopes at two different ages
# (62.15 vs 62.87). Writing the uncentred model as a + b t + g A + d (t A), the
# centred coefficients are beta_1(c) = a + g c and beta_2(c) = b + d c, so gamma
# and delta are invariant to the centring constant while beta_1 and beta_2 move
# with it. The balanced-sample mean is used as the single reference age.
AGE_C <- mean(bal$z_age[bal$wave == 0])
SAMPLES <- list(
  balanced = list(d = bal,  center = AGE_C),
  all      = list(d = alls, center = AGE_C)
)
# descriptive only: the all-available sample's own mean baseline-visit age
AGE_MEAN_ALL <- mean(alls$z_age[!duplicated(alls$id)])

design <- function(d, center, model) {
  t <- d$X_time; Ac <- d$z_age - center
  F <- cbind(beta1 = 1, beta2 = t)
  if (model >= 1) F <- cbind(F, gamma = Ac)
  if (model >= 2) F <- cbind(F, delta = t * Ac)
  F
}

# ---------------------------------------------------------------------------
# 1. Fits
# ---------------------------------------------------------------------------
fits <- list()
for (s in names(SAMPLES)) for (m in 0:2) {
  d <- SAMPLES[[s]]$d
  fits[[paste(s, m)]] <- era_gee_components(d[, Y], design(d, SAMPLES[[s]]$center, m), d$id)
}
main <- fits[["balanced 2"]]

# ---------------------------------------------------------------------------
# 2. Checks
# ---------------------------------------------------------------------------
chk <- list()
ref <- ERA_GEE_fixed(bal[, Y], bal$X_time, rep("gaussian", 4), "independence", y_scale = "raw")
chk$M0_vs_ERA_GEE_fixed <- max(abs(fits[["balanced 0"]]$B - ref$B))
Ac <- bal$z_age - AGE_C; se_diff <- 0
for (q in Y) {
  lf <- lm(bal[[q]] ~ bal$X_time + Ac + I(bal$X_time * Ac))
  s  <- sqrt(diag(vcovCL(lf, cluster = bal$id, type = "HC0", cadjust = FALSE)))
  se_diff <- max(se_diff, abs(s - main$SE[, q]))
}
chk$M2_SE_vs_sandwich <- se_diff
chk$balanced_beta_M2_vs_M0 <- max(abs(main$B[1:2, ] - fits[["balanced 0"]]$B))

# ---------------------------------------------------------------------------
# 3. Main model: age-specific Wave-1 level and per-wave change
#    level(a)  = beta_1q + gamma_q (a - AGE_C)
#    change(a) = beta_2q + delta_q (a - AGE_C)
# ---------------------------------------------------------------------------
lincomb <- function(fit, q, comps, a) {
  nm <- paste0(comps, ":", q)
  c(est = sum(a * fit$B[comps, q]), se = sqrt(as.numeric(t(a) %*% fit$V[nm, nm] %*% a)))
}
OFFS <- c(-10, 0, 10)
agespec <- do.call(rbind, lapply(Y, function(q) do.call(rbind, lapply(OFFS, function(o) {
  l <- lincomb(main, q, c("beta1", "gamma"), c(1, o))
  s <- lincomb(main, q, c("beta2", "delta"), c(1, o))
  data.frame(outcome = YL[[q]], age = AGE_C + o, level = l[1], level_se = l[2],
             change = s[1], change_se = s[2])
}))))
# difference in per-wave change between ages mean+10 and mean-10 = 20 delta
dd20 <- data.frame(outcome = YL[Y], diff = 20 * main$B["delta", Y], se = 20 * main$SE["delta", Y])
dd20$p <- 2 * pnorm(-abs(dd20$diff / dd20$se))

# ---------------------------------------------------------------------------
# 4. Descriptive wave means and the linear summary (balanced)
#    curvature mu0 - 2 mu1 + mu2 is 0 when the two wave-to-wave changes are equal
# ---------------------------------------------------------------------------
Fs  <- cbind(w0 = as.numeric(bal$wave == 0), w1 = as.numeric(bal$wave == 1),
             w2 = as.numeric(bal$wave == 2))
sat <- era_gee_components(bal[, Y], Fs, bal$id)
lin <- do.call(rbind, lapply(Y, function(q) {
  nm <- paste0(c("w0", "w1", "w2"), ":", q)
  mu <- sat$B[, q]; V <- sat$V[nm, nm]
  ctr <- rbind(c(-1, 1, 0), c(0, -1, 1), c(1, -2, 1))
  est <- as.vector(ctr %*% mu); se <- sqrt(diag(ctr %*% V %*% t(ctr)))
  fit0 <- fits[["balanced 0"]]$B[, q]
  data.frame(outcome = YL[[q]], mu0 = mu[1], mu1 = mu[2], mu2 = mu[3],
             fit0 = fit0[1], fit1 = fit0[1] + fit0[2], fit2 = fit0[1] + 2 * fit0[2],
             d01 = est[1], d01_se = se[1], d12 = est[2], d12_se = se[2],
             curv = est[3], curv_se = se[3], curv_p = 2 * pnorm(-abs(est[3] / se[3])))
}))
Cc  <- kronecker(diag(4), t(c(1, -2, 1)))
nmS <- as.vector(sapply(Y, function(q) paste0(c("w0", "w1", "w2"), ":", q)))
cb  <- as.vector(Cc %*% as.vector(sapply(Y, function(q) sat$B[, q])))
Vc  <- Cc %*% sat$V[nmS, nmS] %*% t(Cc)
curv_W <- as.numeric(t(cb) %*% solve(Vc) %*% cb)
curv_p <- pchisq(curv_W, 4, lower.tail = FALSE)

# ---------------------------------------------------------------------------
# 5. Data characteristics that bound the interpretation (balanced persons)
#    Wave 1 PA: baseline "30 minute telephone interview" (startdate_MCQ)
#    Wave 2 PA: FUP1 In-Home Questionnaire (startdate_COF1)
#    Wave 3 PA: FUP2 Questionnaire (startdate_COF2)
#    Age: AGE_NMBR_COM, i.e. at the CLSA baseline visit (startdate_COM)
# ---------------------------------------------------------------------------
rd <- function(f, cols)
  read.csv(file.path(DATA, f), colClasses = "character", check.names = FALSE)[, c("entity_id", cols)]
x0 <- rd("2310007_UofManitoba_SKim_Baseline_CoPv7_Qx_CANUE_PA.csv", c("startdate_COM", "startdate_MCQ"))
x1 <- rd("2310007_UofManitoba_SKim_FUP1_CoPv5_Qx_CANUE_PA.csv", "startdate_COF1")
x2 <- rd("2310007_UofManitoba_SKim_FUP2_CoPv2_Qx_PA.csv", "startdate_COF2")
ids <- unique(bal$entity_id)
dt  <- function(x) as.Date(substr(x, 1, 10))
DC <- dt(x0$startdate_COM[match(ids, x0$entity_id)])
D0 <- dt(x0$startdate_MCQ[match(ids, x0$entity_id)])
D1 <- dt(x1$startdate_COF1[match(ids, x1$entity_id)])
D2 <- dt(x2$startdate_COF2[match(ids, x2$entity_id)])
yrs <- function(a, b) as.numeric(b - a) / 365.25
qs  <- c(.1, .25, .5, .75, .9)
gaps <- rbind(`CLSA baseline visit -> Wave 1 PA` = quantile(yrs(DC, D0), qs, na.rm = TRUE),
              `Wave 1 -> Wave 2` = quantile(yrs(D0, D1), qs, na.rm = TRUE),
              `Wave 2 -> Wave 3` = quantile(yrs(D1, D2), qs, na.rm = TRUE))
yr_rng <- sapply(list(baseline_visit = DC, wave1 = D0, wave2 = D1, wave3 = D2),
                 function(x) paste(range(format(x, "%Y"), na.rm = TRUE), collapse = "-"))
w3_late <- mean(D2 >= as.Date("2020-03-15"), na.rm = TRUE)

# ---------------------------------------------------------------------------
# 6. Report
# ---------------------------------------------------------------------------
con <- file(file.path(OUT, "ERA_GEE_FOUNDATIONAL.md"), open = "wt")
w  <- function(...) writeLines(paste0(...), con)
fp <- function(p) ifelse(p < 1e-4, formatC(p, format = "e", digits = 1), sprintf("%.4f", p))
ci <- function(e, se) sprintf("[%+.5f, %+.5f]", e - 1.96 * se, e + 1.96 * se)
coef_table <- function(f) {
  w("| component | outcome | estimate | SE | 95% CI | p |"); w("|---|---|---|---|---|---|")
  for (k in rownames(f$B)) for (q in Y) {
    e <- f$B[k, q]; se <- f$SE[k, q]
    w("| ", k, " | ", YL[[q]], " | ", sprintf("%+.5f", e), " | ", sprintf("%.5f", se), " | ",
      ci(e, se), " | ", fp(2 * pnorm(-abs(e / se))), " |")
  }
  w("")
  w("Joint Wald tests across the four outcomes (4 df each):")
  w("")
  for (k in rownames(f$B)[-1]) {
    t4 <- era_wald(f, paste0(k, ":", Y))
    w("- ", k, " = 0 for all outcomes: W = ", sprintf("%.1f", t4["W"]), ", p = ", fp(t4["p"]))
  }
  w("")
}

w("# Foundational ERA-GEE model (no partitioning)")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `Analysis/15_era_gee_foundational.R`")
w("")
w("## Definitions")
w("")
w("- **t = 0, 1, 2 denotes the data-collection wave index rather than elapsed calendar time.**")
w("  beta_2q represents the average linear change in Y_q associated with moving from one")
w("  data-collection wave to the next. It is a single linear summary of the overall pattern")
w("  across the three waves, not an annual change and not a constant true slope.")
w("- **A_c = age at the CLSA baseline visit**, centred at the balanced-sample mean of ",
  sprintf("%.2f", AGE_C), " years. beta_1q is the fitted value of Y_q at Wave 1 for a")
w("  participant whose age at the CLSA baseline visit was ", sprintf("%.2f", AGE_C), " years.")
w("- **Main age-adjusted model: M2**, g(mu_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c).")
w("  It evaluates age-related heterogeneity in longitudinal change: delta_q is the difference")
w("  in the per-wave change associated with a one-year difference in baseline-visit age.")
w("  M0 (no age) and M1 (age main effect only) are kept as reference fits; the comparison")
w("  documents why the age x wave term is needed.")
w("- Outcomes are 1-4 frequency categories analysed with a Gaussian family and identity")
w("  link, i.e. as a continuous approximation. Working independence; person-clustered")
w("  sandwich standard errors; all component weights fixed at 1.")
w("")
w("## Samples")
w("")
for (s in names(SAMPLES)) {
  d <- SAMPLES[[s]]$d
  w("- **", s, "**: ", format(length(unique(d$id)), big.mark = ","), " persons, ",
    format(nrow(d), big.mark = ","), " person-waves; age centred at ",
    sprintf("%.2f", SAMPLES[[s]]$center), " (rows per wave: ", paste(table(d$wave), collapse = " / "), ")")
}
w("")
w("## Checks")
w("")
w("| check | max abs difference |"); w("|---|---|")
w("| M0 coefficients vs `ERA_GEE_fixed()` | ", formatC(chk$M0_vs_ERA_GEE_fixed, format = "e", digits = 1), " |")
w("| M2 SEs vs `sandwich::vcovCL` (HC0, no cluster adjustment) | ", formatC(chk$M2_SE_vs_sandwich, format = "e", digits = 1), " |")
w("| balanced: beta1, beta2 in M2 vs M0 | ", formatC(chk$balanced_beta_M2_vs_M0, format = "e", digits = 1), " |")
w("")

w("## Main model M2 · balanced sample")
w("")
coef_table(main)
w("### Model-implied Wave-1 level and per-wave change by baseline-visit age")
w("")
w("| outcome | baseline-visit age | Wave-1 level (SE) | per-wave change (SE) |"); w("|---|---|---|---|")
for (i in seq_len(nrow(agespec))) with(agespec[i, ], w("| ", outcome, " | ", sprintf("%.2f", age), " | ",
  sprintf("%.4f (%.4f)", level, level_se), " | ", sprintf("%+.4f (%.4f)", change, change_se), " |"))
w("")
w("Difference in per-wave change, age ", sprintf("%.2f", AGE_C + 10), " minus age ",
  sprintf("%.2f", AGE_C - 10), " (= 20 x delta):")
w("")
w("| outcome | difference | SE | p |"); w("|---|---|---|---|")
for (i in seq_len(nrow(dd20))) with(dd20[i, ], w("| ", outcome, " | ", sprintf("%+.4f", diff), " | ",
  sprintf("%.4f", se), " | ", fp(p), " |"))
w("")

w("## Reference fits (balanced): why the age x wave term is needed")
w("")
w("### M0 - intercept + wave"); w(""); coef_table(fits[["balanced 0"]])
w("### M1 - + baseline-visit age"); w(""); coef_table(fits[["balanced 1"]])
w("In the balanced sample beta1 and beta2 are identical across M0, M1 and M2: a centred,")
w("time-invariant covariate (and its product with the wave index) is orthogonal to [1, t]")
w("when every person contributes all three waves. gamma in M1 equals the M2 age association")
w("at the middle wave (gamma_M1 = gamma_M2 + delta_M2 x 1).")
w("")

w("## Observed wave means and the linear summary (balanced)")
w("")
w("| outcome | observed W1 / W2 / W3 | M0 fitted W1 / W2 / W3 | change W1->2 (SE) | change W2->3 (SE) | difference of changes (SE) | p |")
w("|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(lin))) with(lin[i, ], w("| ", outcome, " | ",
  sprintf("%.3f / %.3f / %.3f", mu0, mu1, mu2), " | ", sprintf("%.3f / %.3f / %.3f", fit0, fit1, fit2),
  " | ", sprintf("%+.4f (%.4f)", d01, d01_se), " | ", sprintf("%+.4f (%.4f)", d12, d12_se),
  " | ", sprintf("%+.4f (%.4f)", curv, curv_se), " | ", fp(curv_p), " |"))
w("")
w("Joint test that the two wave-to-wave changes are equal for all four outcomes: W = ",
  sprintf("%.1f", curv_W), ", df = 4, p = ", fp(curv_p), ". The equal-change constraint of")
w("the linear specification is therefore a modelling summary, not a description of the data.")
w("")

w("## Data characteristics that bound the interpretation (balanced persons)")
w("")
w("| interval (years) | p10 | p25 | median | p75 | p90 |"); w("|---|---|---|---|---|---|")
for (r in rownames(gaps)) w("| ", r, " | ", paste(sprintf("%.2f", gaps[r, ]), collapse = " | "), " |")
w("")
w("- Calendar years: CLSA baseline visit ", yr_rng["baseline_visit"], "; Wave 1 PA ",
  yr_rng["wave1"], "; Wave 2 ", yr_rng["wave2"], "; Wave 3 ", yr_rng["wave3"], ".")
w("- Wave 1 PA items come from the baseline 30-minute telephone interview (`_MCQ`);")
w("  Wave 2 from the FUP1 In-Home Questionnaire (`_COF1`); Wave 3 from the FUP2")
w("  Questionnaire (`_COF2`).")
w("- ", sprintf("%.1f%%", 100 * w3_late), " of Wave 3 measurements were collected on or after 15 March 2020.")
w("- The balanced sample conditions on completing all three waves (selection/attrition).")
w("")

w("## Sensitivity: all available person-waves")
w("")
w("### M2"); w(""); coef_table(fits[["all 2"]])
w("beta2 across specifications:")
w("")
w("| sample / model | ", paste(YL[Y], collapse = " | "), " |"); w("|---|---|---|---|---|")
for (k in c("balanced 2", "all 0", "all 1", "all 2"))
  w("| ", k, " | ", paste(sprintf("%+.4f", fits[[k]]$B["beta2", Y]), collapse = " | "), " |")
w("")
close(con)

# console summary
cat("checks:\n"); print(unlist(chk))
cat("\nmain model M2 (balanced):\n"); print(round(main$B, 5)); cat("SE:\n"); print(round(main$SE, 5))
for (k in c("gamma", "delta")) { t4 <- era_wald(main, paste0(k, ":", Y)); cat(sprintf("Wald %s: W=%.1f p=%.3g\n", k, t4["W"], t4["p"])) }
cat("\nage-specific level and per-wave change:\n"); print(agespec, digits = 4, row.names = FALSE)
cat("\n20 x delta:\n"); print(dd20, digits = 4, row.names = FALSE)
cat("\nlinear summary:\n"); print(lin[, c("outcome", "mu0", "mu1", "mu2", "fit0", "fit1", "fit2", "d01", "d12", "curv", "curv_p")], digits = 4, row.names = FALSE)
cat("joint equal-change test: W =", round(curv_W, 1), " p =", signif(curv_p, 3), "\n")
cat("\ngaps:\n"); print(round(gaps, 2)); print(yr_rng); cat("W3 on/after 15 Mar 2020:", round(w3_late, 3), "\n")
cat("\nbeta2 sensitivity:\n"); print(round(sapply(c("balanced 2", "all 0", "all 1", "all 2"), function(k) fits[[k]]$B["beta2", ]), 4))
cat("\nall-sample M2:\n"); print(round(fits[["all 2"]]$B, 5)); print(round(fits[["all 2"]]$SE, 5))
cat("\nwritten: DataPrep/out/ERA_GEE_FOUNDATIONAL.md\n")
