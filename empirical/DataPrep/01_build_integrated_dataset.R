# ---------------------------------------------------------------------------
# 01_build_integrated_dataset.R
#
# Build the INTEGRATED longitudinal dataset: every variable that might be
# needed, pulled from the three CLSA Comprehensive files, merged on entity_id,
# and reshaped to person-wave long format.
#
# DELIBERATELY NO CLEANING.
#   Original CLSA codes are preserved exactly as they appear in the source
#   files. Nothing is recoded, nothing is converted to NA, no case is dropped.
#   CLSA uses different missing-value codes in different variables and in
#   different waves, and several of the partitioning variables need an
#   operational definition that has to be decided by looking at the data.
#   Those decisions come AFTER this step, informed by 02_variable_inventory.R.
#
# The only things added are structural: `wave` (0, 1, 2), `X_time` (= wave),
# `X_const` (= 1), and `rowid`.
#
# Output (local only; never committed):
#   DataPrep/out/clsa_integrated_wide.rds
#   DataPrep/out/clsa_integrated_long.rds
#   DataPrep/out/clsa_integrated_long.csv
#
# Run from the project root:
#     Rscript DataPrep/01_build_integrated_dataset.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr)})

DATA <- "CLSA Datafiles for analysis"
OUT  <- "DataPrep/out"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

F_BL   <- file.path(DATA, "2310007_UofManitoba_SKim_Baseline_CoPv7_Qx_CANUE_PA.csv")
F_FUP1 <- file.path(DATA, "2310007_UofManitoba_SKim_FUP1_CoPv5_Qx_CANUE_PA.csv")
F_FUP2 <- file.path(DATA, "2310007_UofManitoba_SKim_FUP2_CoPv2_Qx_PA.csv")

line <- function(ch = "-") cat(strrep(ch, 78), "\n")

# ---------------------------------------------------------------------------
# Variable map.
#
# NOTE the spelling inconsistency across waves, which a naive script gets
# wrong: light-exercise DURATION is PA2_LSPRTHR_MCQ at baseline but
# PA2_LSRTHR_COF1 / PA2_LSRTHR_COF2 at follow-up (LSRTHR, not LSPRTHR).
#
# Physical activity, all three waves:
#   *_freq  participation / frequency item
#   *_hr    hours item, asked conditionally on the frequency item
#
# Time-varying covariates: OWN_OWN, WEA_INCNEEDS, SDC_URBAN_RURAL, AGE_NMBR.
# Baseline-only (time-invariant by design): SEX_ASK_COM, SDC_DRES_COM,
#   SDC_DCGT_COM, ED_UDR04_COM, ED_HIGH_COM, SDC_COB_COM.
# Follow-up sex items (SDC_CURRSEX, SDC_BTHSEX) are carried through so that
#   consistency with the baseline item can be checked later.
# ---------------------------------------------------------------------------

VARS <- list(
  bl = c(
    entity_id       = "entity_id",
    walk_freq       = "PA2_WALK_MCQ",     walk_hr    = "PA2_WALKHR_MCQ",
    lsport_freq     = "PA2_LSPRT_MCQ",    lsport_hr  = "PA2_LSPRTHR_MCQ",
    msport_freq     = "PA2_MSPRT_MCQ",    msport_hr  = "PA2_MSPRTHR_MCQ",
    ssport_freq     = "PA2_SSPRT_MCQ",    ssport_hr  = "PA2_SSPRTHR_MCQ",
    own             = "OWN_OWN_COM",
    incneeds        = "WEA_INCNEEDS_MCQ",
    urban_rural     = "SDC_URBAN_RURAL_COM",
    age             = "AGE_NMBR_COM",
    # baseline-only, time-invariant
    sex_ask         = "SEX_ASK_COM",
    yrs_in_canada   = "SDC_DRES_COM",
    cultural_bg     = "SDC_DCGT_COM",
    educ4           = "ED_UDR04_COM",
    educ_high       = "ED_HIGH_COM",
    country_birth   = "SDC_COB_COM",
    age_grp         = "AGE_GRP_COM"
  ),
  f1 = c(
    entity_id       = "entity_id",
    walk_freq       = "PA2_WALK_COF1",    walk_hr    = "PA2_WALKHR_COF1",
    lsport_freq     = "PA2_LSPRT_COF1",   lsport_hr  = "PA2_LSRTHR_COF1",
    msport_freq     = "PA2_MSPRT_COF1",   msport_hr  = "PA2_MSPRTHR_COF1",
    ssport_freq     = "PA2_SSPRT_COF1",   ssport_hr  = "PA2_SSPRTHR_COF1",
    own             = "OWN_OWN_COF1",
    incneeds        = "WEA_INCNEEDS_COF1",
    urban_rural     = "SDC_URBAN_RURAL_COF1",
    age             = "AGE_NMBR_COF1",
    sex_current     = "SDC_CURRSEX_COF1",
    sex_birth       = "SDC_BTHSEX_COF1"
  ),
  f2 = c(
    entity_id       = "entity_id",
    walk_freq       = "PA2_WALK_COF2",    walk_hr    = "PA2_WALKHR_COF2",
    lsport_freq     = "PA2_LSPRT_COF2",   lsport_hr  = "PA2_LSRTHR_COF2",
    msport_freq     = "PA2_MSPRT_COF2",   msport_hr  = "PA2_MSPRTHR_COF2",
    ssport_freq     = "PA2_SSPRT_COF2",   ssport_hr  = "PA2_SSPRTHR_COF2",
    own             = "OWN_OWN_COF2",
    incneeds        = "WEA_INCNEEDS_COF2",
    urban_rural     = "SDC_URBAN_RURAL_COF2",
    age             = "AGE_NMBR_COF2",
    sex_current     = "SDC_CURRSEX_COF2"
  )
)

# every value is read as character so that no CLSA code is silently coerced,
# reformatted, or turned into NA by the reader
read_wave <- function(path, map, label) {
  hdr <- names(read_csv(path, n_max = 0, show_col_types = FALSE,
                        name_repair = "minimal"))
  missing_cols <- setdiff(unname(map), hdr)
  if (length(missing_cols))
    stop(label, ": columns not found in the source file: ",
         paste(missing_cols, collapse = ", "))

  d <- read_csv(path, col_select = all_of(unname(map)),
                col_types = cols(.default = col_character()),
                show_col_types = FALSE, progress = FALSE)
  d <- d[, unname(map), drop = FALSE]
  names(d) <- names(map)
  cat(sprintf("  %-8s %6d rows x %2d vars   (%s)\n",
              label, nrow(d), ncol(d), basename(path)))
  d
}

line("="); cat("READING SOURCE FILES (all values kept as character)\n"); line("=")
bl <- read_wave(F_BL,   VARS$bl, "baseline")
f1 <- read_wave(F_FUP1, VARS$f1, "FUP1")
f2 <- read_wave(F_FUP2, VARS$f2, "FUP2")

# ---------------------------------------------------------------------------
# Wide merge. Full join so that nobody is lost: a person present in only one
# wave is retained, with NA for the waves in which they do not appear.
# Those NAs are STRUCTURAL (person not observed), and are distinct from the
# CLSA missing codes carried inside the variables themselves.
# ---------------------------------------------------------------------------
line("="); cat("MERGING ON entity_id (full join, nobody dropped)\n"); line("=")

suffixed <- function(d, sfx) {
  id <- d[["entity_id"]]
  d2 <- d[, setdiff(names(d), "entity_id"), drop = FALSE]
  names(d2) <- paste0(names(d2), sfx)
  cbind(entity_id = id, d2, stringsAsFactors = FALSE)
}

wide <- suffixed(bl, "_w0") |>
  full_join(suffixed(f1, "_w1"), by = "entity_id") |>
  full_join(suffixed(f2, "_w2"), by = "entity_id")

cat(sprintf("  distinct entity_id: baseline %d | FUP1 %d | FUP2 %d | merged %d\n",
            n_distinct(bl$entity_id), n_distinct(f1$entity_id),
            n_distinct(f2$entity_id), n_distinct(wide$entity_id)))
cat(sprintf("  wide dataset: %d rows x %d columns\n", nrow(wide), ncol(wide)))

cat("\n  wave presence pattern (person appears in a wave's source file):\n")
pres <- data.frame(
  w0 = !is.na(wide$walk_freq_w0) | !is.na(wide$sex_ask_w0),
  w1 = wide$entity_id %in% f1$entity_id,
  w2 = wide$entity_id %in% f2$entity_id)
pres$w0 <- wide$entity_id %in% bl$entity_id
tb <- table(paste0(as.integer(pres$w0), as.integer(pres$w1), as.integer(pres$w2)))
for (k in names(sort(tb, decreasing = TRUE)))
  cat(sprintf("    waves present %s : %6d persons (%.1f%%)\n",
              k, tb[[k]], 100 * tb[[k]] / nrow(wide)))

# ---------------------------------------------------------------------------
# Long reshape: one row per person-wave.
# ---------------------------------------------------------------------------
line("="); cat("RESHAPING TO PERSON-WAVE LONG\n"); line("=")

TV <- c("walk_freq","walk_hr","lsport_freq","lsport_hr",
        "msport_freq","msport_hr","ssport_freq","ssport_hr",
        "own","incneeds","urban_rural","age","sex_current","sex_birth")
TI <- c("sex_ask","yrs_in_canada","cultural_bg","educ4","educ_high",
        "country_birth","age_grp")

get_col <- function(nm, w) if (nm %in% names(wide)) wide[[nm]] else NA_character_

long <- do.call(rbind, lapply(0:2, function(w) {
  d <- data.frame(entity_id = wide$entity_id, wave = w, stringsAsFactors = FALSE)
  for (v in TV) d[[v]] <- get_col(paste0(v, "_w", w), w)
  for (v in TI) d[[v]] <- wide[[paste0(v, "_w0")]]   # baseline-only, repeated
  d
}))
long <- long[order(long$entity_id, long$wave), ]
long$X_const <- 1
long$X_time  <- long$wave
long$rowid   <- seq_len(nrow(long))
rownames(long) <- NULL

cat(sprintf("  long dataset: %d rows (= %d persons x 3 waves) x %d columns\n",
            nrow(long), n_distinct(long$entity_id), ncol(long)))
cat("  columns:\n    ", paste(names(long), collapse = ", "), "\n", sep = "")

cat("\n  NOTE: rows exist for all three waves for every person. A row whose\n")
cat("  variables are all NA means the person was not observed in that wave;\n")
cat("  a row carrying CLSA codes such as 8, 9, -88888 or -99999 means the\n")
cat("  person WAS observed and the item itself is coded. Those two kinds of\n")
cat("  'missing' are deliberately left distinguishable at this stage.\n")

# ---------------------------------------------------------------------------
saveRDS(wide, file.path(OUT, "clsa_integrated_wide.rds"))
saveRDS(long, file.path(OUT, "clsa_integrated_long.rds"))
write_csv(long, file.path(OUT, "clsa_integrated_long.csv"), na = "")

line("="); cat("WRITTEN\n"); line("=")
for (f in c("clsa_integrated_wide.rds", "clsa_integrated_long.rds",
            "clsa_integrated_long.csv"))
  cat(sprintf("  %-32s %8.1f MB\n", f,
              file.size(file.path(OUT, f)) / 1024^2))
cat("\n  These files contain CLSA data and must stay local. Never commit them.\n")
