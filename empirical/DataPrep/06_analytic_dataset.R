# ---------------------------------------------------------------------------
# 06_analytic_dataset.R
#
# Build the analytic dataset from the integrated + derived data.
#
# TWO RULES APPLIED HERE, both settled from descriptive evidence:
#
# (1) PARTITIONING VARIABLES ARE BASELINE-ONLY AND TIME-INVARIANT.
#     The MOB-ERA-GEE specification requires Z to be person-level. A
#     time-varying Z would let the same person fall on both sides of a split,
#     turning the tree from "which kinds of PEOPLE have different trajectories"
#     into "which PERSON-WAVES show different parameter processes", which is not
#     the subgroup interpretation the paper makes. Baseline values are therefore
#     carried across all waves as `z_*`. Later-wave versions stay in the file
#     for descriptive and sensitivity use.
#
#     Age in particular MUST be baseline age: wave-specific age increases almost
#     deterministically with time, so using it as Z would confound the
#     longitudinal time effect with the subgrouping structure.
#
# (2) A PERSON-WAVE IS RETAINED ONLY IF ALL FOUR OUTCOMES ARE OBSERVED.
#     The implementation treats the Q = 4 outcome vector jointly within a
#     person-wave. Partial outcome missingness is rare (146 person-waves, 0.18%
#     of observed rows, 143 persons) and arises only after CLSA-declared item
#     non-response codes are set to NA; the PA block itself was administered
#     all-or-nothing.
#
# Outputs (local only):
#   DataPrep/out/analytic_long_all.rds       every retained person-wave
#   DataPrep/out/analytic_long_balanced.rds  persons complete at all three waves
#
# Run from the project root:
#     Rscript DataPrep/06_analytic_dataset.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

OUT <- "DataPrep/out"
d   <- readRDS(file.path(OUT, "clsa_integrated_long_v2.rds"))
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
line <- function(ch = "-") cat(strrep(ch, 78), "\n")

line("="); cat("BUILDING THE ANALYTIC DATASET\n"); line("=")

# ---------------------------------------------------------------------------
# 1. Baseline-only partitioning variables, carried across waves
# ---------------------------------------------------------------------------
bl <- d[d$wave == 0, ]
CUTS <- c(10, 15, 20, 25, 30, 35, 40)

zb <- data.frame(entity_id = bl$entity_id, stringsAsFactors = FALSE)
zb$z_sex        <- bl$sex
zb$z_racialized <- bl$racialized
zb$z_eth8       <- bl$eth8
zb$z_educ4      <- bl$educ4f
zb$z_educ3      <- bl$educ3f
zb$z_own        <- bl$own3          # baseline dwelling ownership
zb$z_incneeds   <- bl$incneeds5     # baseline income adequacy
zb$z_urban      <- bl$urban2        # baseline urban/rural (code 9 -> Urban)
zb$z_urban3     <- bl$urban3
zb$z_age        <- bl$age_n         # BASELINE age
zb$z_immigrant  <- factor(ifelse(is.na(bl$is_immigrant), NA,
                          ifelse(bl$is_immigrant == 1, "Immigrant", "Non-immigrant")),
                          levels = c("Non-immigrant", "Immigrant"))
zb$z_yrs_imm    <- bl$yrs_imm       # NA for non-immigrants: immigrants-only use
for (k in CUTS) {
  zb[[paste0("z_imm_cut", k)]] <- bl[[paste0("imm_status_cut", k)]]
  zb[[paste0("z_g4_cut", k)]]  <- bl[[paste0("g4_cut", k)]]
}
cat(sprintf("  baseline Z variables: %d\n", ncol(zb) - 1))
cat("    z_sex, z_racialized, z_eth8, z_educ4/3, z_own, z_incneeds, z_urban(3),\n")
cat("    z_age, z_immigrant, z_yrs_imm, z_imm_cut{...}, z_g4_cut{...}\n")

d <- d[, setdiff(names(d), names(zb)[-1])]
d <- merge(d, zb, by = "entity_id", all.x = TRUE, sort = FALSE)
d <- d[order(d$entity_id, d$wave), ]

# ---------------------------------------------------------------------------
# 2. Retain a person-wave only if all four outcomes are observed
# ---------------------------------------------------------------------------
n_obs_rows <- sum(rowSums(!is.na(d[, Y])) > 0 | !is.na(d$own) | !is.na(d$age))
keep <- rowSums(!is.na(d[, Y])) == 4L
cat(sprintf("\n  person-wave rows before: %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("  retained (all 4 outcomes observed): %s\n", format(sum(keep), big.mark = ",")))
cat(sprintf("  dropped, partial outcomes (1-3 of 4): %s\n",
            format(sum(rowSums(!is.na(d[, Y])) %in% 1:3), big.mark = ",")))
cat(sprintf("  dropped, no outcomes / person not in wave: %s\n",
            format(sum(rowSums(!is.na(d[, Y])) == 0), big.mark = ",")))
a <- d[keep, ]

# ---------------------------------------------------------------------------
# 3. Completeness flags and the balanced subsample
# ---------------------------------------------------------------------------
nw <- table(a$entity_id)
a$n_waves   <- as.integer(nw[as.character(a$entity_id)])
a$complete3 <- a$n_waves == 3L
a$rowid     <- seq_len(nrow(a))          # required by the mob fit factories
a$id        <- a$entity_id               # clustering variable

bal <- a[a$complete3, ]
bal$rowid <- seq_len(nrow(bal))

pw <- tapply(a$wave, a$entity_id, function(x) paste(sort(x), collapse = ""))
cat(sprintf("\n  persons with >=1 retained wave: %s\n",
            format(length(unique(a$entity_id)), big.mark = ",")))
cat("  wave pattern among retained person-waves:\n")
tb <- sort(table(pw), decreasing = TRUE)
for (k in names(tb))
  cat(sprintf("    waves {%s}: %6s persons (%.1f%%)\n", k,
              format(tb[[k]], big.mark = ","), 100 * tb[[k]] / length(pw)))

cat(sprintf("\n  ALL (unbalanced): %s rows, %s persons\n",
            format(nrow(a), big.mark = ","), format(length(unique(a$entity_id)), big.mark = ",")))
cat(sprintf("  BALANCED (3 complete waves): %s rows, %s persons\n",
            format(nrow(bal), big.mark = ","), format(length(unique(bal$entity_id)), big.mark = ",")))

# ---------------------------------------------------------------------------
saveRDS(a,   file.path(OUT, "analytic_long_all.rds"))
saveRDS(bal, file.path(OUT, "analytic_long_balanced.rds"))
line("="); cat("WRITTEN\n"); line("=")
cat("  analytic_long_all.rds       (every retained person-wave)\n")
cat("  analytic_long_balanced.rds  (persons complete at all three waves)\n")
cat("  CLSA data -- stays local, never committed\n")
