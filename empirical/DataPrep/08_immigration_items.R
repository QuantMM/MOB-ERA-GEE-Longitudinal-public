# ---------------------------------------------------------------------------
# 08_immigration_items.R
#
# Stage: RECORD + INVENTORY for the immigration-defining items.
# Nothing here is a cleaning or coding decision.
#
# STARTING POINT. Everyone in the integrated data (every baseline participant),
# with their outcome availability, and then the baseline items that identify an
# immigrant. Values are recorded raw (character) and only their RESPONSE STATUS
# is classified; what people actually answered (which country, how many years,
# what age) is left for a later stage.
#
# THE ITEMS
#   SDC_COB_COM      "In what country were you born?"            ASKED
#                    001 = Canada; 998 / 999 / 777 = CLSA missing
#   SDC_YACA_YR_COM  "In what year did you first come to Canada   ASKED
#                    to live?"   9998 / 9999 / 7777 = CLSA missing
#   SDC_DRES_COM     Length of time in Canada since immigration   DERIVED
#   SDC_DAIM_COM     Age at time of immigration                   DERIVED
#                    both: 996 = excluded, non-immigrant; 999 = not answered
#   SDC_FIMM_COM     CLSA Immigration Flag (1 / 2 / 9)            DERIVED
#                    recorded as a cross-check only
#
#   Only COB and YACA are questions put to the participant; DRES and DAIM are
#   CLSA derivations from them. All of these exist at baseline only (absent from
#   the FUP1 and FUP2 files), so they are person-level by construction.
#
# WHAT "INDICATES IMMIGRANT" MEANS HERE
#   COB   a valid country other than Canada
#   DRES  a valid number of years (not 996, not 999)
#   DAIM  a valid age (not 996, not 999)
#   Each item is classified separately; the report shows how the three agree.
#
# DISCLOSURE. Cells < 6 are suppressed. This is an internal working document:
# some suppressed cells can be recovered from the margins, so complementary
# suppression is required before any of these tables leave the project.
#
# Run from the project root:
#     Rscript DataPrep/08_immigration_items.R
# Writes (local only):
#     DataPrep/out/immigration_items_raw.rds
#     DataPrep/out/IMMIGRATION_ITEMS.md
# ---------------------------------------------------------------------------

OUT  <- "DataPrep/out"
DATA <- "CLSA Datafiles for analysis"
BASE <- file.path(DATA, "2310007_UofManitoba_SKim_Baseline_CoPv7_Qx_CANUE_PA.csv")
ITEMS <- c("SDC_COB_COM", "SDC_YACA_YR_COM", "SDC_DRES_COM", "SDC_DAIM_COM",
           "SDC_FIMM_COM", "SDC_GCB_COM")
Y <- c("y_walk", "y_lsport", "y_msport", "y_ssport")

sup <- function(n) ifelse(n > 0 & n < 6, "<6", format(n, big.mark = ",", trim = TRUE))
pct <- function(n, d) ifelse(n > 0 & n < 6, "", sprintf("%.1f%%", 100 * n / d))

# ---------------------------------------------------------------------------
# 1. Everyone, with outcome availability by wave
# ---------------------------------------------------------------------------
L <- readRDS(file.path(OUT, "clsa_integrated_long_v2.rds"))
L$n_y <- rowSums(!is.na(L[, Y]))
per <- data.frame(entity_id = sort(unique(L$entity_id)), stringsAsFactors = FALSE)
for (w in 0:2) {
  s <- L[L$wave == w, c("entity_id", "n_y")]
  per[[paste0("ny_w", w)]] <- s$n_y[match(per$entity_id, s$entity_id)]
}
cw <- as.matrix(per[, paste0("ny_w", 0:2)])
per$complete_waves <- rowSums(!is.na(cw) & cw == 4L)
per$any_outcome    <- rowSums(!is.na(cw) & cw > 0L) > 0

# ---------------------------------------------------------------------------
# 2. The immigration items, raw
# ---------------------------------------------------------------------------
hdr  <- strsplit(gsub('"', "", readLines(BASE, n = 1)), ",")[[1]]
keep <- c("entity_id", ITEMS)
stopifnot(all(keep %in% hdr))
raw <- read.csv(BASE, colClasses = "character", check.names = FALSE)[, keep]
raw[] <- lapply(raw, trimws)
saveRDS(raw, file.path(OUT, "immigration_items_raw.rds"))

d <- merge(per, raw, by = "entity_id", all.x = TRUE, sort = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))

cob  <- num(d$SDC_COB_COM); yaca <- num(d$SDC_YACA_YR_COM)
dres <- num(d$SDC_DRES_COM); daim <- num(d$SDC_DAIM_COM)

d$st_cob <- ifelse(d$SDC_COB_COM == "" | is.na(d$SDC_COB_COM), "blank",
            ifelse(cob %in% c(998, 999, 777), paste0("missing ", cob),
            ifelse(cob == 1, "Canada", "other country")))
d$st_yaca <- ifelse(d$SDC_YACA_YR_COM == "" | is.na(d$SDC_YACA_YR_COM), "blank",
             ifelse(yaca %in% c(9998, 9999, 7777), paste0("missing ", yaca), "valid year"))
d$st_dres <- ifelse(dres == 996, "996 non-immigrant", ifelse(dres == 999, "999 not answered", "valid years"))
d$st_daim <- ifelse(daim == 996, "996 non-immigrant", ifelse(daim == 999, "999 not answered", "valid age"))
d$st_fimm <- c(`1` = "1 immigrant", `2` = "2 not an immigrant", `9` = "9 not answered")[d$SDC_FIMM_COM]

d$ind_cob  <- d$st_cob  == "other country"
d$ind_dres <- d$st_dres == "valid years"
d$ind_daim <- d$st_daim == "valid age"
d$n_ind    <- d$ind_cob + d$ind_dres + d$ind_daim
d$imm_any  <- d$n_ind >= 1
d$imm_all  <- d$n_ind == 3

# ---------------------------------------------------------------------------
# 3. Report
# ---------------------------------------------------------------------------
con <- file(file.path(OUT, "IMMIGRATION_ITEMS.md"), open = "wt")
w <- function(...) writeLines(paste0(...), con)
tab <- function(x, lab, N = nrow(d)) {
  t <- table(x, useNA = "ifany"); names(t)[is.na(names(t))] <- "(NA)"
  w("| ", lab, " | persons | % |"); w("|---|---|---|")
  for (k in names(t)) w("| ", k, " | ", sup(t[[k]]), " | ", pct(t[[k]], N), " |")
  w("")
}
N <- nrow(d)

w("# Immigration-defining items — response inventory")
w("")
w("**Generated:** ", format(Sys.Date(), "%d %B %Y"), " · `DataPrep/08_immigration_items.R`")
w("")
w("Stage: record and inventory. No coding decision is taken here. Cells < 6 are")
w("suppressed; this is an internal document and needs complementary suppression")
w("before any table leaves the project.")
w("")
w("## 1. Starting population: everyone, with outcome availability")
w("")
w("All **", format(N, big.mark = ","), "** baseline participants in the integrated data.")
w("Waves with all four physical-activity outcomes observed:")
w("")
tab(factor(d$complete_waves, levels = 0:3, labels = paste(0:3, "complete waves")), "outcome availability")

w("## 2. The items")
w("")
w("| item | source | CLSA label |")
w("|---|---|---|")
w("| `SDC_COB_COM` | **asked**: \"In what country were you born?\" | Country of birth |")
w("| `SDC_YACA_YR_COM` | **asked**: \"In what year did you first come to Canada to live?\" | Year arrival in Canada |")
w("| `SDC_DRES_COM` | derived by CLSA | Length of time in Canada since immigration |")
w("| `SDC_DAIM_COM` | derived by CLSA | Age at Time of Immigration |")
w("| `SDC_FIMM_COM` | derived by CLSA (cross-check only) | Immigration Flag |")
w("")
w("All five are baseline-only; none appears in the FUP1 or FUP2 files.")
w("")
w("### Response status, each item separately")
w("")
tab(d$st_cob, "`SDC_COB_COM`"); tab(d$st_yaca, "`SDC_YACA_YR_COM`")
tab(d$st_dres, "`SDC_DRES_COM`"); tab(d$st_daim, "`SDC_DAIM_COM`")
tab(d$st_fimm, "`SDC_FIMM_COM`")

w("## 3. Agreement among the three immigrant indications")
w("")
w("Indication = COB a country other than Canada; DRES a valid number of years;")
w("DAIM a valid age.")
w("")
pat <- paste0(ifelse(d$ind_cob, "COB", "—"), " · ", ifelse(d$ind_dres, "DRES", "—"),
              " · ", ifelse(d$ind_daim, "DAIM", "—"))
tp <- sort(table(pat), decreasing = TRUE)
w("| pattern (COB · DRES · DAIM) | persons |"); w("|---|---|")
for (k in names(tp)) w("| ", k, " | ", sup(tp[[k]]), " |")
w("")
w("- indicated by **at least one** item: **", sup(sum(d$imm_any)), "**")
w("- indicated by **all three**: **", sup(sum(d$imm_all)), "**")
w("- indicated by one or two but not all three: **", sup(sum(d$imm_any & !d$imm_all)), "**")
w("")
w("### Where the partial cases sit on the asked questions")
w("")
pc <- d[d$imm_any & !d$imm_all, ]
if (nrow(pc)) {
  x <- table(COB = pc$st_cob, YACA = pc$st_yaca)
  w("| COB \\ YACA | ", paste(colnames(x), collapse = " | "), " |")
  w("|", paste(rep("---", ncol(x) + 1), collapse = "|"), "|")
  for (i in rownames(x)) w("| ", i, " | ", paste(sup(x[i, ]), collapse = " | "), " |")
  w("")
}

w("## 4. Cross-checks")
w("")
w("Against the CLSA Immigration Flag:")
w("")
x <- table(flag = d$st_fimm, indicated = ifelse(d$imm_all, "all three",
           ifelse(d$imm_any, "one or two", "none")))
w("| FIMM \\ indications | ", paste(colnames(x), collapse = " | "), " |")
w("|", paste(rep("---", ncol(x) + 1), collapse = "|"), "|")
for (i in rownames(x)) w("| ", i, " | ", paste(sup(x[i, ]), collapse = " | "), " |")
w("")
w("Against the pipeline's current `is_immigrant` (from `SDC_DRES_COM` alone):")
w("")
cur <- L[L$wave == 0, c("entity_id", "is_immigrant")]
d$cur <- cur$is_immigrant[match(d$entity_id, cur$entity_id)]
x <- table(current = ifelse(is.na(d$cur), "NA", ifelse(d$cur == 1, "immigrant", "non-immigrant")),
           any_indication = ifelse(d$imm_any, "immigrant", "not indicated"))
w("| current \\ any indication | ", paste(colnames(x), collapse = " | "), " |")
w("|", paste(rep("---", ncol(x) + 1), collapse = "|"), "|")
for (i in rownames(x)) w("| ", i, " | ", paste(sup(x[i, ]), collapse = " | "), " |")
w("")

w("## 5. The immigrant group against outcome availability")
w("")
w("| | everyone | any outcome | ≥1 complete wave | 3 complete waves (balanced) |")
w("|---|---|---|---|---|")
row <- function(lab, s) w("| ", lab, " | ", sup(sum(s)), " | ", sup(sum(s & d$any_outcome)), " | ",
                          sup(sum(s & d$complete_waves >= 1)), " | ", sup(sum(s & d$complete_waves == 3)), " |")
row("all participants", rep(TRUE, N))
row("immigrant — any indication", d$imm_any)
row("immigrant — all three", d$imm_all)
row("not indicated as immigrant", !d$imm_any)
w("")
w("Share of the immigrant group (any indication) retained at each step:")
w("")
a <- sum(d$imm_any); b <- sum(!d$imm_any)
w("| | immigrant | not indicated |"); w("|---|---|---|")
w("| ≥1 complete wave | ", sprintf("%.1f%%", 100 * sum(d$imm_any & d$complete_waves >= 1) / a), " | ",
  sprintf("%.1f%%", 100 * sum(!d$imm_any & d$complete_waves >= 1) / b), " |")
w("| 3 complete waves | ", sprintf("%.1f%%", 100 * sum(d$imm_any & d$complete_waves == 3) / a), " | ",
  sprintf("%.1f%%", 100 * sum(!d$imm_any & d$complete_waves == 3) / b), " |")
w("")
close(con)

cat("written: DataPrep/out/IMMIGRATION_ITEMS.md\n")
cat("written: DataPrep/out/immigration_items_raw.rds (raw, character, local only)\n")
