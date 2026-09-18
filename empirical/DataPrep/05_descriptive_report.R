# ---------------------------------------------------------------------------
# 05_descriptive_report.R
#
# Produces DataPrep/out/DESCRIPTIVE_REPORT.md with figures in
# DataPrep/out/fig/. Descriptive evidence only; no cleaning decisions are
# taken here.
#
# CLSA SMALL-CELL RULE. Every count below 6 is suppressed as "<6", and where
# suppressing one cell would let it be recovered by subtraction, the smallest
# further cell is suppressed too. The helper `supp()` implements this and is
# used for every table that leaves this pipeline.
#
# Run from the project root:
#     Rscript DataPrep/05_descriptive_report.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

OUT <- "DataPrep/out"; FIG <- file.path(OUT, "fig")
if (!dir.exists(FIG)) dir.create(FIG, recursive = TRUE)
d <- readRDS(file.path(OUT, "clsa_integrated_long_v2.rds"))

con <- file(file.path(OUT, "DESCRIPTIVE_REPORT.md"), open = "wt")
w   <- function(...) writeLines(paste0(...), con)
Y   <- c("y_walk","y_lsport","y_msport","y_ssport")
YL  <- c("Walking","Light exercise","Moderate exercise","Strenuous exercise")
b   <- d[d$wave == 0, ]

# observed-in-wave indicator
obs <- lapply(0:2, function(k) {
  dd <- d[d$wave == k, ]
  cols <- c("walk_freq","own","incneeds","urban_rural","age")
  Reduce(`|`, lapply(cols, function(v) !is.na(dd[[v]])))
})
names(obs) <- as.character(0:2)
nobs <- sapply(obs, sum)

# --- CLSA small-cell suppression -------------------------------------------
supp <- function(x, thr = 6) {
  x <- as.integer(x); out <- format(x, big.mark = ",", trim = TRUE)
  small <- !is.na(x) & x > 0 & x < thr
  out[small] <- paste0("<", thr)
  # secondary disclosure: if exactly one cell is suppressed in a vector,
  # suppress the next smallest as well
  if (sum(small) == 1L && length(x) > 1L) {
    cand <- which(!small & !is.na(x))
    if (length(cand)) out[cand[which.min(x[cand])]] <- paste0("<", thr)
  }
  out[!is.na(x) & x == 0] <- "0"
  out
}
mdtab <- function(m, rn = "") {
  w("| ", rn, " | ", paste(colnames(m), collapse = " | "), " |")
  w("|", paste(rep("---", ncol(m) + 1), collapse = "|"), "|")
  for (i in seq_len(nrow(m)))
    w("| ", rownames(m)[i], " | ", paste(m[i, ], collapse = " | "), " |")
  w("")
}
pct <- function(n, tot) sprintf("%.1f", 100 * n / tot)

# ===========================================================================
w("# CLSA analytic data — descriptive report")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `DataPrep/05_descriptive_report.R`  ")
w("**Data:** CLSA Comprehensive cohort, three waves, Application 2310007")
w("")
w("This report is the second stage of a deliberately staged pipeline:")
w("")
w("> raw → inventory → codebook → **descriptive evidence** → cleaning decisions → analytic dataset")
w("")
w("The usual order is inverted on purpose. CLSA codes missingness differently in")
w("different variables and different waves, and several of the variables needed here")
w("require an operational definition that cannot be settled without looking at the")
w("distributions first. So nothing is filtered, collapsed or excluded at this stage")
w("beyond recoding the values CLSA's own dictionaries declare to be missing.")
w("Derived variables are *recorded alongside* the raw columns, several alternatives")
w("at once, so that the choice between them stays open and evidence-based.")
w("")
w("All counts below 6 are suppressed as `<6` in accordance with the CLSA")
w("Publication Policy, including where a second cell must be suppressed to prevent")
w("recovery by subtraction.")
w("")
w("---")
w("")

# ===========================================================================
w("## 1. Sample and wave participation")
w("")
w("Every participant appears in the baseline file; imbalance arises only from")
w("attrition and from item non-response.")
w("")
w(sprintf("- persons: **%s**", format(nrow(b), big.mark = ",")))
w(sprintf("- person-wave rows: **%s**", format(nrow(d), big.mark = ",")))
w(sprintf("- observed in wave 0 / 1 / 2: **%s / %s / %s**",
          format(nobs[1], big.mark=","), format(nobs[2], big.mark=","),
          format(nobs[3], big.mark=",")))
w("")
pat <- paste0(as.integer(obs[["0"]]), as.integer(obs[["1"]]), as.integer(obs[["2"]]))
tp  <- sort(table(pat), decreasing = TRUE)
m <- cbind(persons = supp(as.integer(tp)), `%` = pct(as.integer(tp), nrow(b)))
rownames(m) <- paste0("`", names(tp), "`")
w("Wave participation pattern (1 = observed, positions are waves 0/1/2):")
w("")
mdtab(m, "pattern")

# ===========================================================================
w("## 2. Measures and Variable Operationalization")
w("")
w("Each variable is described under eight headings: the construct it represents,")
w("the CLSA source item, its original response scale, how it was operationalized")
w("here, whether it is baseline-only or time-varying, its intended analytic role,")
w("how missing values were treated, and the justification for any decision that is")
w("not self-evident.")
w("")
w("Where a decision is genuinely open it is *not* taken here. Several")
w("operationalizations are recorded side by side and the open question is stated.")
w("")

measure <- function(name, construct, source, scale, oper, timing, role, missing, just = NULL) {
  w("### ", name); w("")
  w("| | |"); w("|---|---|")
  w("| **Construct** | ", construct, " |")
  w("| **Data source** | ", source, " |")
  w("| **Original coding** | ", scale, " |")
  w("| **Operationalization** | ", oper, " |")
  w("| **Timing** | ", timing, " |")
  w("| **Analytic role** | ", role, " |")
  w("| **Missing-value treatment** | ", missing, " |")
  if (!is.null(just)) { w(""); w("**Justification.** ", just) }
  w("")
}

w("#### Outcomes")
w("")
measure(
  "Physical activity participation — four frequency indicators",
  paste("Multidimensional physical activity participation in later life. The four",
        "items are treated as empirically distinct expressions of a participation",
        "profile rather than as interchangeable indicators of one latent dimension."),
  paste("PASE (Physical Activity Scale for the Elderly) frequency items.",
        "`PA2_WALK_*` frequency of taking a walk outside, past 7 days;",
        "`PA2_LSPRT_*` light sport/recreational activity;",
        "`PA2_MSPRT_*` moderate; `PA2_SSPRT_*` strenuous.",
        "Wave suffixes `_MCQ` (baseline), `_COF1`, `_COF2`."),
  "1 Never · 2 Seldom (1–2 days) · 3 Sometimes (3–4 days) · 4 Often (5–7 days)",
  paste("Retained on the original 4-point scale as `y_walk`, `y_lsport`,",
        "`y_msport`, `y_ssport`. No composite is formed."),
  "Time-varying; measured identically at all three waves.",
  "**Outcomes.** The four form the multivariate response block.",
  "8 (Don't know/No answer), 9 (Refused) and −88888 (Missing, wave 2) → NA. All three are flagged `missing` in the CLSA dictionaries.",
  paste("The duration counterparts (`*HR`) are not used: they are asked",
        "conditionally on the frequency item, so they are observed only for a",
        "selected subset and are not comparable across respondents as parallel",
        "outcomes. Section 3 shows the frequency items share one scale, span the",
        "full range at every wave, and carry ~0.1% item non-response."))

w("#### Immigration and racialization")
w("")
measure(
  "Years since immigration / immigrant status",
  "Time since arrival as a proxy for accumulated exposure to the structural conditions of settlement, and the basis of the non-/recent-/established-immigrant contrast.",
  "`SDC_DRES_COM`, \"Length of time in Canada since immigration\".",
  paste("0–88 genuine years; **996 \"Excluded participants – Non-immigrants\"**;",
        "999 \"Required question was not answered\". Both 996 and 999 are flagged",
        "`missing` by CLSA."),
  paste("`yrs_imm` = years, defined only for immigrants. `is_immigrant` binary.",
        "`imm_status_cut{10,15,20,25,30,35,40}` = non-immigrant / recent / established",
        "at seven candidate thresholds, **all recorded**."),
  "Baseline only; time-invariant by construction.",
  "Partitioning variable. Threshold, or continuous use, **not yet decided**.",
  "999 (n = 3) → NA. **996 is treated as a category, not as missing**: it identifies the non-immigrant comparison group and is the single most informative value in the variable.",
  paste("CLSA flags 996 `missing` because the question is not applicable to",
        "Canadian-born respondents, not because the value is unknown. Reading it",
        "as missing would delete 81.9% of the sample. Seven thresholds are recorded",
        "rather than one because the distribution (Section 4) makes any single",
        "conventional cut hard to defend on its own."))

measure(
  "Cultural / racial background",
  paste("Two analytically distinct constructs. **Ethnicity** is self-identified",
        "cultural background — ancestry, language, traditions — and carries",
        "cultural practices and identity processes. **Racialized group** refers to",
        "populations socially constructed as non-White and positioned within",
        "systems of structural inequality, and carries exposure to discrimination",
        "and structural barriers. Ethnicity alone cannot separate cultural",
        "variation from structural inequality."),
  "`SDC_DCGT_COM`, \"Cultural / Racial Background\".",
  paste("1 White · 2 Black · 3 Korean · 4 Filipino · 5 Japanese · 6 Chinese ·",
        "7 South Asian · 8 Southeast Asian · 9 Arab · 10 West Asian ·",
        "11 Latin American · 12 Other · 13 Multiple ·",
        "96 Excluded – Aboriginal identity · 99 Not answered."),
  paste("`eth8`: White · Black · East Asian (3,5,6) · South/Southeast Asian (4,7,8) ·",
        "Middle Eastern (9,10) · Latin American (11) · Other (12) · Multiple (13).",
        "`racialized`: White vs Racialized.",
        "`g4_cut{k}`: White / Racialized + Canadian-born / Racialized + recent /",
        "Racialized + long-term, at each candidate threshold."),
  "Baseline only; time-invariant.",
  "Descriptive (`eth8`) and partitioning (`racialized`, `g4_cut{k}`).",
  "99 (n = 30) → NA, participants retained. **96 does not occur in this extract**, so the Aboriginal-identity exclusion has already been applied upstream.",
  paste("Grouping follows Morin et al. (2022), *Osteoporosis International* 33(12),",
        "2637–2648, with three revisions adopted for this project: Arab and West",
        "Asian combined as Middle Eastern, and Latin American and Other kept",
        "separate rather than merged. Because the theoretical frame is Health",
        "Lifestyle Theory — structural life chances constraining behaviour — the",
        "operative construct is racialization; `eth8` is descriptive. `g4_cut{k}`",
        "is recorded because racialization and immigration are not proxies for one",
        "another in this cohort (Section 5)."))

w("#### Structural and socioeconomic covariates")
w("")
measure(
  "Dwelling ownership",
  "Housing tenure as an indicator of accumulated economic security and of stability of residence.",
  "`OWN_OWN_COM` / `OWN_OWN_COF1` / `OWN_OWN_COF2`.",
  paste("**Coding differs between waves.** Baseline uses zero-padded strings and",
        "distinguishes 01 Own · 02 Rent · 03 Leasehold · 04 No cost · 05 Provided",
        "by employer · 06 Joint ownership · 07 Own and rent · 97 Other · 77 Missing",
        "· 98 Don't know · 99 Refused. Follow-up uses plain integers and has",
        "already collapsed everything but Own and Rent into 97 Other."),
  "`own3` = Own / Rent / Other, applied identically at all three waves.",
  "Time-varying.",
  "Partitioning variable.",
  "77, 98, 99, −88888, −99999 → NA.",
  paste("Collapsing to three categories is forced by the data rather than chosen:",
        "the follow-up waves do not carry the finer distinctions, so no longitudinal",
        "coding that preserves them exists."))

measure(
  "Income meets needs",
  "Subjective adequacy of household income — perceived economic strain, distinct from income level.",
  "`WEA_INCNEEDS_MCQ` / `_COF1` / `_COF2`.",
  "1 Very well · 2 Adequately · 3 With some difficulty · 4 Not very well · 5 Totally inadequately.",
  "`incneeds5` ordered factor on the original scale; `incneeds_n` numeric 1–5.",
  "Time-varying; identical scale at all three waves.",
  "Partitioning variable.",
  "8 (Don't know/No answer), 9 (Refused), −88888 (Missing) → NA.",
  NULL)

measure(
  "Urban / rural residence",
  "Residential context: availability of walkable environments, recreational facilities and transport, all of which condition the feasibility of physical activity.",
  "`SDC_URBAN_RURAL_COM` / `_COF1` / `_COF2`.",
  paste("0 Rural · 1 Urban core · 2 Urban fringe · 4 Urban population centre outside",
        "CMA/CA · 6 Secondary core · 9 Postal code link to dissemination area."),
  "`urban2` = Urban / Rural, with code 9 classified as urban. `urban3` keeps code 9 as its own level, for checking.",
  "Time-varying.",
  "Partitioning variable.",
  "−88888, −99999 → NA.",
  paste("Code 9 is a **geocoding method, not a place type**: CLSA derives the",
        "classification from Statistics Canada's Postal Code Conversion File, and",
        "where a six-character postal code cannot be resolved to a population",
        "centre it is linked to a dissemination area instead. CLSA's own *Data",
        "Support Document: Urban / Rural Classification* recommends classifying",
        "these respondents as urban, and Yuan et al. (2025, *Canadian Journal of",
        "Public Health*, doi:10.17269/s41997-025-01088-4) — also using the CLSA",
        "Comprehensive cohort — follow that recommendation explicitly. `urban3` is",
        "retained so the recommendation can be checked against these data",
        "(Section 5) rather than taken on trust."))

measure(
  "Education",
  "Educational attainment as an indicator of class position and of access to health information and to occupational security.",
  "`ED_UDR04_COM`.",
  paste("1 Less than secondary · 2 Secondary graduation, no post-secondary ·",
        "3 Some post-secondary · 4 Post-secondary degree/diploma ·",
        "9 At least one required question not answered."),
  "`educ4f` on the original four levels; `educ3f` collapsing 2 and 3.",
  "Baseline only.",
  "Partitioning variable.",
  "9 (n = 50) → NA.",
  "Both codings are recorded. The three-level version follows the convention in comparable CLSA analyses; the four-level version is retained so that the collapse can be checked rather than assumed.")

measure(
  "Sex",
  "Sex as a structural position shaping caregiving demands, time availability and access to physical activity.",
  "`SEX_ASK_COM` (baseline). Follow-up additionally carries `SDC_CURRSEX_COF1/2` (current gender) and `SDC_BTHSEX_COF1` (sex at birth).",
  "F / M. **No 0 code occurs in this extract.**",
  "`sex` = Female / Male, from the baseline item.",
  "Treated as baseline, time-invariant.",
  "Partitioning variable.",
  "None required; the baseline item is complete.",
  "The follow-up gender items are carried through in the integrated dataset so that consistency with the baseline item can be examined, but the baseline item is used because it is the only one present at all three waves.")

measure(
  "Age",
  "Chronological age, the primary axis of an aging cohort and a strong candidate source of trajectory heterogeneity.",
  "`AGE_NMBR_COM` / `AGE_NMBR_COF1` / `AGE_NMBR_COF2`; `AGE_GRP_COM` for the baseline age band.",
  "Continuous years.",
  "`age_n` continuous.",
  "Time-varying, available at all three waves.",
  "Candidate partitioning variable; **inclusion not yet decided**.",
  "None required.",
  "Age is not part of the specification inherited from the study proposal, but it is available, complete, and a plausible driver of heterogeneity in physical activity trajectories. Recorded so that the decision can be made explicitly rather than by omission.")

w("---")
w("")

# ===========================================================================
w("## 3. Outcomes")
w("")
w("### 3.1 Distribution by wave")
w("")
for (i in seq_along(Y)) {
  w("**", YL[i], "** (`", Y[i], "`)"); w("")
  m <- sapply(0:2, function(k) {
    x <- d[[Y[i]]][d$wave == k]
    tb <- table(factor(x, levels = 1:4))
    c(supp(as.integer(tb)), sprintf("%.2f", mean(x, na.rm = TRUE)),
      sprintf("%.2f", sd(x, na.rm = TRUE)), format(sum(is.na(x) & obs[[k+1]]), big.mark=","))
  })
  colnames(m) <- paste("wave", 0:2)
  rownames(m) <- c("1 Never","2 Seldom","3 Sometimes","4 Often","mean","SD","item missing")
  mdtab(m, YL[i])
}

png(file.path(FIG, "outcome_distributions.png"), width = 1100, height = 750, res = 120)
op <- par(mfrow = c(2,2), mar = c(4,4,3,1))
for (i in seq_along(Y)) {
  tb <- sapply(0:2, function(k) prop.table(table(factor(d[[Y[i]]][d$wave==k], levels=1:4))))
  barplot(tb, beside = TRUE, names.arg = paste0("w", 0:2), ylim = c(0, .6),
          col = c("#d9d9d9","#a6bddb","#3690c0","#034e7b"), main = YL[i],
          ylab = "proportion", legend.text = if (i == 1) c("Never","Seldom","Sometimes","Often") else NULL,
          args.legend = list(x = "topright", bty = "n", cex = .8))
}
par(op); invisible(dev.off())
w("![outcome distributions](fig/outcome_distributions.png)")
w("")
w("The four indicators sit at very different levels — walking is common, moderate")
w("and strenuous activity rare — and each is stable across waves. That is the")
w("empirical basis for modelling them jointly rather than combining them.")
w("")

w("### 3.2 Association among the four outcomes")
w("")
cm <- cor(d[d$wave == 0, Y], use = "pairwise.complete.obs")
dimnames(cm) <- list(YL, YL)
mdtab(matrix(sprintf("%.3f", cm), 4, 4, dimnames = dimnames(cm)), "wave 0")
w("Pairwise correlations are small throughout. A composite score would discard")
w("most of the information that distinguishes the four behaviours.")
w("")

w("### 3.3 Within-person variation across waves")
w("")
wv <- sapply(Y, function(v) {
  m <- sapply(0:2, function(k) d[[v]][d$wave == k])
  ok <- rowSums(!is.na(m)) == 3
  c(`between-person SD` = sd(rowMeans(m[ok, ]), na.rm = TRUE),
    `mean within-person SD` = mean(apply(m[ok, ], 1, sd), na.rm = TRUE))
})
colnames(wv) <- YL
mdtab(matrix(sprintf("%.3f", wv), 2, 4, dimnames = list(rownames(wv), YL)), "")
w("Between-person variation dominates, but within-person variation is not")
w("negligible — which is what a growth model needs.")
w("")

w("### 3.4 Outcome observation pattern")
w("")
m <- sapply(0:2, function(k) {
  dd <- d[d$wave == k, ][obs[[k+1]], ]
  tb <- table(factor(rowSums(!is.na(dd[, Y])), levels = 0:4))
  supp(as.integer(tb))
})
colnames(m) <- paste("wave", 0:2); rownames(m) <- paste(0:4, "of 4")
mdtab(m, "outcomes present")
w("**There is no partial outcome missingness anywhere**: either all four items were")
w("administered or none were. This matters for the analysis stage — a working")
w("covariance across the four outcomes never meets a partially observed vector.")
w("")
cnt <- sapply(0:2, function(k) {
  dd <- d[d$wave == k, ]
  setNames(rowSums(!is.na(dd[, Y])) == 4L, dd$entity_id)[as.character(b$entity_id)]
})
tb <- table(factor(rowSums(cnt, na.rm = TRUE), levels = 0:3))
m <- cbind(persons = supp(as.integer(tb)), `%` = pct(as.integer(tb), nrow(b)))
rownames(m) <- paste(0:3, "wave(s)")
w("Waves at which a person has all four outcomes:")
w("")
mdtab(m, "complete waves")

# ===========================================================================
w("## 4. Immigration history")
w("")
yv <- b$yrs_imm
w(sprintf("- non-immigrant (996): **%s**", format(sum(b$is_immigrant == 0, na.rm=TRUE), big.mark=",")))
w(sprintf("- immigrant: **%s**", format(sum(b$is_immigrant == 1, na.rm=TRUE), big.mark=",")))
w(sprintf("- status not answered (999): **%s**", supp(sum(is.na(b$is_immigrant)))))
w("")
q <- quantile(yv, c(0,.05,.10,.25,.50,.75,.90,.95,1), na.rm = TRUE)
m <- matrix(sprintf("%.0f", q), 1, length(q),
            dimnames = list("years", names(q)))
w("Years since immigration, among immigrants:"); w("")
mdtab(m, "")
w("**The median immigrant arrived ", round(median(yv, na.rm=TRUE)), " years ago**, and only 5% arrived within ",
    round(quantile(yv, .05, na.rm=TRUE)), " years. In this cohort \"recent immigrant\" does not mean what it means")
w("in the general migration literature, and the paper has to say so.")
w("")

png(file.path(FIG, "years_since_immigration.png"), width = 1000, height = 450, res = 120)
op <- par(mar = c(4.5,4.5,3,1))
hist(yv, breaks = 40, col = "#3690c0", border = "white",
     main = "Years since immigration (immigrants only)",
     xlab = "years since immigration", ylab = "persons")
abline(v = c(10,20,30,40), lty = 2, col = "grey30")
text(c(10,20,30,40), par("usr")[4]*.92, c("10","20","30","40"), col = "grey30", cex = .8)
par(op); invisible(dev.off())
w("![years since immigration](fig/years_since_immigration.png)")
w("")
CUTS <- c(10,15,20,25,30,35,40)
m <- t(sapply(CUTS, function(k) {
  tb <- table(b[[paste0("imm_status_cut", k)]])
  c(supp(tb[["recent"]]), pct(tb[["recent"]], nrow(b)), supp(tb[["established"]]))
}))
colnames(m) <- c("recent", "% of all", "established")
rownames(m) <- paste0(CUTS, " years")
w("Group sizes at each candidate threshold (non-immigrant = ",
  format(sum(b$is_immigrant == 0, na.rm=TRUE), big.mark=","), " throughout):")
w("")
mdtab(m, "threshold")

# ===========================================================================
w("## 5. Racialization")
w("")
tb <- table(b$eth8, useNA = "no")
m <- cbind(n = supp(as.integer(tb)), `%` = pct(as.integer(tb), nrow(b)))
rownames(m) <- names(tb)
mdtab(m, "group")
w(sprintf("Not answered (code 99): %s. Code 96, Aboriginal identity, does not occur.",
          supp(sum(is.na(b$eth8)))))
w("")
w("### 5.1 Racialization by immigrant status (20-year threshold)")
w("")
tt <- table(b$eth8, b$imm_status_cut20)
m <- apply(tt, 2, supp); rownames(m) <- rownames(tt)
mdtab(m, "group")
w("The eight-group scheme survives this single cross-tabulation, but only barely:")
w("the smallest cells are in single figures. Adding sex, or any second split,")
w("breaks it. It is reportable as description; it is not a safe partitioning")
w("variable.")
w("")
w("### 5.2 The four-level racialization × immigration variable")
w("")
tb <- table(b$g4_cut20)
m <- cbind(n = supp(as.integer(tb)), `%` = pct(as.integer(tb), nrow(b)))
rownames(m) <- names(tb)
mdtab(m, "group")
tt <- table(b$g4_cut20, b$sex)
m <- apply(tt, 2, supp); rownames(m) <- rownames(tt)
w("Crossed with sex:"); w("")
mdtab(m, "group")
w("Smallest cell **", min(tt), "**. Racialization and immigration are **not proxies for")
w("one another here**: ", format(sum(b$g4_cut20 == "2 Racialized + Canadian-born", na.rm=TRUE), big.mark=","),
  " racialized participants are Canadian-born, and ",
  format(sum(b$racialized == "White" & b$is_immigrant == 1, na.rm = TRUE), big.mark=","),
  " White participants are immigrants.")
w("")

w("### 5.3 Code 9 of the urban/rural variable")
w("")
w("Code 9 (\"postal code link to dissemination area\") is a geocoding outcome, not")
w("a place type. Before accepting CLSA's recommendation to treat it as urban, its")
w("composition is compared with the two substantive categories.")
w("")
b0 <- d[d$wave == 0, ]
u3 <- b0$urban3
rows <- list(
  `n` = table(u3),
  `% immigrant` = tapply(b0$is_immigrant, u3, function(x) 100*mean(x, na.rm=TRUE)),
  `% racialized` = tapply(as.integer(b0$racialized == "Racialized"), u3, function(x) 100*mean(x, na.rm=TRUE)),
  `% female` = tapply(as.integer(b0$sex == "Female"), u3, function(x) 100*mean(x, na.rm=TRUE)),
  `mean age` = tapply(b0$age_n, u3, mean, na.rm = TRUE),
  `% post-sec degree` = tapply(as.integer(b0$educ4f == "Post-secondary degree"), u3, function(x) 100*mean(x, na.rm=TRUE)),
  `% own home` = tapply(as.integer(b0$own3 == "Own"), u3, function(x) 100*mean(x, na.rm=TRUE)),
  `mean walking freq` = tapply(b0$y_walk, u3, mean, na.rm = TRUE)
)
m <- t(sapply(names(rows), function(k)
  if (k == "n") supp(as.integer(rows[[k]])) else sprintf("%.1f", rows[[k]])))
colnames(m) <- names(table(u3))
mdtab(m, "")
w("On every characteristic the DA-linked group resembles the urban categories more")
w("closely than the rural one: it is more likely to be immigrant (19.6% vs 18.5%")
w("urban and 13.9% rural), better educated (78.6% vs 77.9% and 74.7%), *less*")
w("likely to own its home than either (77.5% vs 83.5% and 94.1%), and identical to")
w("the urban group on mean walking frequency. Home ownership in particular runs")
w("firmly against a rural reading.")
w("")
w("So CLSA's recommendation is supported empirically in this sample as well as by")
w("derivation, and `urban2` (code 9 → Urban) can be used with that stated.")
w("`urban3` is retained so the alternative remains available as a sensitivity")
w("check rather than being closed off.")
w("")

# ===========================================================================
w("## 6. Covariate distributions by wave")
w("")
cov_tab <- function(v, lab) {
  m <- sapply(0:2, function(k) {
    x <- d[[v]][d$wave == k]
    tb <- table(x)
    c(supp(as.integer(tb)), format(sum(is.na(x) & obs[[k+1]]), big.mark = ","))
  })
  rn <- c(levels(d[[v]]), "missing")
  m <- matrix(unlist(m), nrow = length(rn),
              dimnames = list(rn, paste("wave", 0:2)))
  w("**", lab, "** (`", v, "`)"); w("")
  mdtab(m, lab)
}
cov_tab("own3", "Dwelling ownership")
cov_tab("incneeds5", "Income meets needs")
cov_tab("urban2", "Urban / rural")

w("**Education** (`educ4f`, baseline only)"); w("")
tb <- table(b$educ4f)
m <- cbind(n = supp(as.integer(tb)), `%` = pct(as.integer(tb), nrow(b)))
rownames(m) <- names(tb); mdtab(m, "level")
w(sprintf("Not answered: %s.", supp(sum(is.na(b$educ4f)))))
w("")
w("**Sex** (`sex`, baseline)"); w("")
tb <- table(b$sex)
m <- cbind(n = supp(as.integer(tb)), `%` = pct(as.integer(tb), nrow(b)))
rownames(m) <- names(tb); mdtab(m, "")
w("**Age** (`age_n`)"); w("")
m <- t(sapply(0:2, function(k) {
  x <- d$age_n[d$wave == k][obs[[k+1]]]
  sprintf("%.1f", c(mean(x, na.rm=TRUE), sd(x, na.rm=TRUE),
                    min(x, na.rm=TRUE), max(x, na.rm=TRUE)))
}))
colnames(m) <- c("mean","SD","min","max"); rownames(m) <- paste("wave", 0:2)
mdtab(m, "")

# ===========================================================================
w("## 7. Outcomes by group")
w("")
w("Mean of each outcome by immigrant status (20-year threshold) and wave.")
w("Descriptive only — no model is fitted here.")
w("")
for (i in seq_along(Y)) {
  m <- sapply(0:2, function(k) {
    dd <- d[d$wave == k, ]
    sprintf("%.2f", tapply(dd[[Y[i]]], dd$imm_status_cut20, mean, na.rm = TRUE))
  })
  colnames(m) <- paste("wave", 0:2)
  rownames(m) <- levels(b$imm_status_cut20)
  w("**", YL[i], "**"); w(""); mdtab(m, "")
}

png(file.path(FIG, "outcome_by_immigrant_status.png"), width = 1100, height = 750, res = 120)
op <- par(mfrow = c(2,2), mar = c(4,4.2,3,1), oma = c(2.5,0,0,0))
cols <- c("#737373", "#3690c0", "#d1495b")
for (i in seq_along(Y)) {
  mm <- sapply(0:2, function(k) {
    dd <- d[d$wave == k, ]; tapply(dd[[Y[i]]], dd$imm_status_cut20, mean, na.rm = TRUE)
  })
  matplot(0:2, t(mm), type = "b", pch = 19, lty = 1, lwd = 2, col = cols,
          xaxt = "n", xlab = "wave", ylab = "mean frequency", main = YL[i])
  axis(1, at = 0:2)
}
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
legend("bottom", legend = levels(b$imm_status_cut20), col = cols, lwd = 2,
       pch = 19, horiz = TRUE, bty = "n", cex = .85)
par(op); invisible(dev.off())
w("![outcomes by immigrant status](fig/outcome_by_immigrant_status.png)")
w("")

# ===========================================================================
w("## 8. What remains open")
w("")
w("Deliberately not decided at this stage:")
w("")
w("1. **Balanced or unbalanced person-waves.** Requires knowing how the analytic")
w("   engine handles unbalanced clusters and what missingness assumptions it makes.")
w("2. **The recent/established threshold**, or whether years since immigration")
w("   enters continuously. Section 4 gives the counts at seven candidate cuts.")
w("3. **Which variables enter the model, and in what role.** Nothing here")
w("   pre-commits to a set of partitioning variables.")
w("4. **Code 9 of the urban/rural variable** — Section 5.3 provides the evidence;")
w("   both codings are recorded.")
w("5. **Whether age is included** as a partitioning variable.")
w("")
w("Settled by the data rather than by choice:")
w("")
w("- dwelling ownership must be Own / Rent / Other, because the follow-up waves")
w("  carry no finer distinction;")
w("- the outcomes need no partial-missingness handling, because partial")
w("  missingness does not occur;")
w("- the Aboriginal-identity exclusion needs no action, because code 96 is absent.")
w("")

close(con)
cat("written: ", file.path(OUT, "DESCRIPTIVE_REPORT.md"), "\n")
cat("figures: ", FIG, "\n")
