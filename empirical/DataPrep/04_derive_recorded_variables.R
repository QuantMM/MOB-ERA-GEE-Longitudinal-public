# ---------------------------------------------------------------------------
# 04_derive_recorded_variables.R
#
# Add derived variables to the integrated dataset. RAW COLUMNS ARE PRESERVED
# UNCHANGED; everything here is appended alongside them.
#
# These are RECORDED, not DECIDED. Creating `imm_status_cut20` does not commit
# the analysis to a 20-year threshold, and creating `g4` does not commit it to
# analysing within those four groups. Several alternatives are recorded side by
# side precisely so that the choice can be made later, from descriptive
# evidence, rather than now.
#
# The one class of recoding done here is NOT a researcher judgement: CLSA's own
# dictionaries carry a `missing` flag, and codes CLSA itself declares missing
# are set to NA in the derived (not the raw) columns. Without that, no
# descriptive statistic can be computed at all. Every such code is listed in
# the report.
#
# Run from the project root:
#     Rscript DataPrep/04_derive_recorded_variables.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

OUT  <- "DataPrep/out"
long <- readRDS(file.path(OUT, "clsa_integrated_long.rds"))

num <- function(x) suppressWarnings(as.numeric(x))
line <- function(ch = "-") cat(strrep(ch, 78), "\n")

line("="); cat("DERIVING RECORDED VARIABLES (raw columns untouched)\n"); line("=")

# ---------------------------------------------------------------------------
# 1. Physical activity outcomes
#    Same 4-point scale in all three waves:
#      1 Never | 2 Seldom (1-2 days) | 3 Sometimes (3-4 days) | 4 Often (5-7 days)
#    CLSA-declared missing: 8 (Don't know/No answer), 9 (Refused),
#    -88888 (Missing, wave 2 only).
# ---------------------------------------------------------------------------
PA_MISS <- c(8, 9, -88888, -99999)
for (v in c("walk", "lsport", "msport", "ssport")) {
  src <- paste0(v, "_freq")
  x <- num(long[[src]])
  x[x %in% PA_MISS] <- NA
  long[[paste0("y_", v)]] <- x
}
cat(sprintf("  y_walk / y_lsport / y_msport / y_ssport : 4-point, %s set to NA\n",
            paste(PA_MISS, collapse = "/")))

# ---------------------------------------------------------------------------
# 2. Years since immigration, and immigrant status at several thresholds
#    996 = "Excluded participants - Non-immigrants"  (a CATEGORY for us, not
#          a missing value: it identifies the comparison group)
#    999 = "Required question was not answered"      (genuinely unknown)
# ---------------------------------------------------------------------------
yr <- num(long$yrs_in_canada)
long$yrs_imm      <- ifelse(yr < 996, yr, NA)          # years, immigrants only
long$is_immigrant <- ifelse(yr == 996, 0L, ifelse(yr == 999, NA, 1L))

CUTS <- c(10, 15, 20, 25, 30, 35, 40)
for (k in CUTS)
  long[[paste0("imm_status_cut", k)]] <- factor(
    ifelse(yr == 996, "non-immigrant",
    ifelse(yr == 999, NA,
    ifelse(yr <  k,   "recent", "established"))),
    levels = c("non-immigrant", "established", "recent"))
cat(sprintf("  yrs_imm, is_immigrant, imm_status_cut{%s}\n",
            paste(CUTS, collapse = ",")))

# ---------------------------------------------------------------------------
# 3. Cultural / racial background
#    Eight descriptive groups, following Morin et al. (2022) with the
#    revisions recorded for this project (Middle Eastern = Arab + West Asian;
#    Latin American and Other kept separate; Multiple kept separate).
#    99 = not answered -> NA.  96 (Aboriginal) is absent from this extract.
# ---------------------------------------------------------------------------
E8 <- c("1"="White", "2"="Black",
        "3"="East Asian", "5"="East Asian", "6"="East Asian",
        "4"="South/Southeast Asian", "7"="South/Southeast Asian",
        "8"="South/Southeast Asian",
        "9"="Middle Eastern", "10"="Middle Eastern",
        "11"="Latin American", "12"="Other", "13"="Multiple")
long$eth8 <- factor(unname(E8[long$cultural_bg]),
                    levels = c("White","Black","East Asian","South/Southeast Asian",
                               "Middle Eastern","Latin American","Other","Multiple"))
long$racialized <- factor(ifelse(is.na(long$eth8), NA,
                                 ifelse(long$eth8 == "White", "White", "Racialized")),
                          levels = c("White", "Racialized"))
cat("  eth8 (8 descriptive groups), racialized (binary)\n")

# ---------------------------------------------------------------------------
# 4. Racialization crossed with immigration history.
#    Recorded at every candidate threshold so the threshold stays open.
# ---------------------------------------------------------------------------
for (k in CUTS) {
  st <- as.character(long[[paste0("imm_status_cut", k)]])
  long[[paste0("g4_cut", k)]] <- factor(
    ifelse(is.na(long$racialized) | is.na(st), NA,
    ifelse(long$racialized == "White", "1 White",
    ifelse(st == "non-immigrant", "2 Racialized + Canadian-born",
    ifelse(st == "recent",        "3 Racialized + recent immigrant",
                                  "4 Racialized + long-term immigrant")))),
    levels = c("1 White", "2 Racialized + Canadian-born",
               "3 Racialized + recent immigrant", "4 Racialized + long-term immigrant"))
}
cat(sprintf("  g4_cut{%s} (racialization x immigration)\n", paste(CUTS, collapse = ",")))

# ---------------------------------------------------------------------------
# 5. Dwelling ownership.
#    Baseline codes are zero-padded strings and distinguish 01-07; follow-up
#    uses plain integers and has already collapsed everything but Own/Rent into
#    97 "Other". A common longitudinal coding therefore has to be Own/Rent/Other.
#    CLSA-declared missing: 98, 99, -99999. Baseline 77 ("Missing") also.
# ---------------------------------------------------------------------------
own_raw <- trimws(long$own)
own_n   <- num(own_raw)
long$own3 <- factor(
  ifelse(is.na(own_n) | own_n %in% c(77, 98, 99, -88888, -99999), NA,
  ifelse(own_n == 1, "Own",
  ifelse(own_n == 2, "Rent", "Other"))),
  levels = c("Own", "Rent", "Other"))
cat("  own3 (Own/Rent/Other; forced by the follow-up coding)\n")

# ---------------------------------------------------------------------------
# 6. Income meets needs. Same 1-5 scale in all waves.
#    CLSA-declared missing: 8, 9, -88888.
# ---------------------------------------------------------------------------
inc <- num(long$incneeds)
inc[inc %in% c(8, 9, -88888, -99999)] <- NA
long$incneeds5 <- factor(inc, levels = 1:5,
  labels = c("Very well","Adequately","With some difficulty",
             "Not very well","Totally inadequately"))
long$incneeds_n <- inc
cat("  incneeds5 (ordered factor), incneeds_n (numeric 1-5)\n")

# ---------------------------------------------------------------------------
# 7. Urban / rural.
#    0 Rural | 1 Urban core | 2 Urban fringe | 4 Urban population centre
#    outside CMA/CA | 6 Secondary core | 9 Postal code link to dissemination area
#
#    Code 9 is a GEOCODING method, not a place type. CLSA's own Data Support
#    Document recommends classifying it as urban; Yuan et al. (2025, CJPH,
#    doi:10.17269/s41997-025-01088-4) follow that recommendation explicitly.
#    Both codings are recorded so the choice can be checked empirically.
# ---------------------------------------------------------------------------
ur <- num(long$urban_rural)
ur[ur %in% c(-88888, -99999)] <- NA
long$urban2 <- factor(ifelse(is.na(ur), NA,
                      ifelse(ur == 0, "Rural", "Urban")),
                      levels = c("Urban", "Rural"))              # 9 -> Urban
long$urban3 <- factor(ifelse(is.na(ur), NA,
                      ifelse(ur == 0, "Rural",
                      ifelse(ur == 9, "DA-linked", "Urban"))),
                      levels = c("Urban", "Rural", "DA-linked")) # 9 kept apart
cat("  urban2 (code 9 -> Urban, per CLSA recommendation), urban3 (9 kept separate)\n")

# ---------------------------------------------------------------------------
# 8. Education, sex, age.
#    educ4: 1 < secondary | 2 secondary, no post-sec | 3 some post-sec
#           4 post-secondary degree/diploma | 9 not answered  [CLSA: MISSING]
# ---------------------------------------------------------------------------
ed <- num(long$educ4); ed[ed == 9] <- NA
long$educ4f <- factor(ed, levels = 1:4,
  labels = c("< secondary","Secondary, no post-sec","Some post-secondary",
             "Post-secondary degree"))
long$educ3f <- factor(ifelse(is.na(ed), NA, ifelse(ed == 1, 1, ifelse(ed == 4, 3, 2))),
  levels = 1:3,
  labels = c("< secondary","Secondary or some post-secondary","Post-secondary degree"))
long$sex <- factor(long$sex_ask, levels = c("F","M"), labels = c("Female","Male"))
long$age_n <- num(long$age)
cat("  educ4f, educ3f, sex, age_n\n")

# ---------------------------------------------------------------------------
saveRDS(long, file.path(OUT, "clsa_integrated_long_v2.rds"))
line("="); cat("WRITTEN\n"); line("=")
cat(sprintf("  clsa_integrated_long_v2.rds   %d rows x %d columns  (%.1f MB)\n",
            nrow(long), ncol(long),
            file.size(file.path(OUT, "clsa_integrated_long_v2.rds")) / 1024^2))
cat("  raw columns preserved; derived columns appended\n")
cat("  CLSA data -- stays local, never committed\n")
