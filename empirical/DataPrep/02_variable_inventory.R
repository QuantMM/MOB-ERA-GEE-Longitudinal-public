# ---------------------------------------------------------------------------
# 02_variable_inventory.R
#
# Variable-by-variable inventory of the integrated dataset, for deciding
# operational definitions. Nothing is recoded here either.
#
# For each variable it reports, separately for each wave:
#   - how many persons were observed in that wave at all
#   - the complete value distribution, every distinct code, with counts
#   - which values look like CLSA reserved / missing codes
#   - structural NA (person not in that wave) vs item-level codes
#
# CLSA uses different reserved codes in different variables and waves
# (8, 9, 77, 96, 97, 98, 99, -88888, -99999, ...), so the point of this
# report is to see what each variable actually contains before deciding.
#
# Run from the project root:
#     Rscript DataPrep/02_variable_inventory.R
# Writes DataPrep/out/variable_inventory.txt as well as printing.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

OUT  <- "DataPrep/out"
long <- readRDS(file.path(OUT, "clsa_integrated_long.rds"))
wide <- readRDS(file.path(OUT, "clsa_integrated_wide.rds"))

con <- file(file.path(OUT, "variable_inventory.txt"), open = "wt")
say <- function(...) { txt <- paste0(...); cat(txt); writeLines(sub("\n$", "", txt), con) }
line <- function(ch = "-") say(strrep(ch, 78), "\n")

# codes CLSA commonly reserves; used only to FLAG values for attention
RESERVED <- c("7","8","9","77","88","96","97","98","99",
              "777","888","999","-8","-88888","-99999",
              "88888","99999","-7777","7777")

# which persons appear in each wave's source file at all
in_wave <- list(
  `0` = !is.na(wide$sex_ask_w0) | !is.na(wide$age_w0) | !is.na(wide$walk_freq_w0),
  `1` = Reduce(`|`, lapply(grep("_w1$", names(wide), value = TRUE),
                           function(v) !is.na(wide[[v]]))),
  `2` = Reduce(`|`, lapply(grep("_w2$", names(wide), value = TRUE),
                           function(v) !is.na(wide[[v]])))
)
n_obs <- sapply(in_wave, sum)

TV <- c("walk_freq","walk_hr","lsport_freq","lsport_hr",
        "msport_freq","msport_hr","ssport_freq","ssport_hr",
        "own","incneeds","urban_rural","age","sex_current","sex_birth")
TI <- c("sex_ask","yrs_in_canada","cultural_bg","educ4","educ_high",
        "country_birth","age_grp")

say("\n"); line("=")
say("CLSA INTEGRATED DATASET -- VARIABLE INVENTORY\n")
line("=")
say(sprintf("\npersons: %d | person-wave rows: %d\n", nrow(wide), nrow(long)))
say(sprintf("observed in wave 0 / 1 / 2: %d / %d / %d\n\n",
            n_obs[["0"]], n_obs[["1"]], n_obs[["2"]]))

# ---------------------------------------------------------------------------
show_var <- function(v, waves, invariant = FALSE) {
  line("="); say(sprintf("%s%s\n", v, if (invariant) "   [baseline only, repeated across waves]" else ""))
  line("=")
  for (w in waves) {
    d  <- long[long$wave == w, ]
    x  <- d[[v]]
    seen <- in_wave[[as.character(w)]][match(d$entity_id, wide$entity_id)]

    n_struct <- sum(is.na(x) &  !seen)     # person not observed in this wave
    n_itemna <- sum(is.na(x) &   seen)     # observed, but this item is blank
    vals <- x[!is.na(x)]
    tb   <- sort(table(vals), decreasing = TRUE)

    say(sprintf("\n  wave %d   observed persons = %d\n", w, sum(seen)))
    say(sprintf("    NA, person not in this wave : %6d\n", n_struct))
    say(sprintf("    NA, observed but item blank : %6d\n", n_itemna))
    say(sprintf("    non-missing values          : %6d   distinct = %d\n",
                length(vals), length(tb)))
    if (!length(tb)) next

    numeric_like <- suppressWarnings(!any(is.na(as.numeric(names(tb)))))
    if (length(tb) > 25 && numeric_like) {
      q <- suppressWarnings(quantile(as.numeric(vals), c(0,.01,.25,.5,.75,.99,1), na.rm = TRUE))
      say(sprintf("    continuous-looking; min/1%%/25%%/50%%/75%%/99%%/max = %s\n",
                  paste(format(q, trim = TRUE), collapse = " / ")))
      flagged <- intersect(names(tb), RESERVED)
      big <- names(tb)[suppressWarnings(abs(as.numeric(names(tb))) >= 900)]
      flagged <- union(flagged, big)
      if (length(flagged)) {
        say("    values needing an explicit decision:\n")
        for (k in flagged)
          say(sprintf("      %-10s %6d  (%.2f%% of observed)\n",
                      k, tb[[k]], 100 * tb[[k]] / sum(seen)))
      }
      say("    ten most frequent: ")
      say(paste(sprintf("%s=%d", names(tb)[1:min(10,length(tb))],
                        as.integer(tb)[1:min(10,length(tb))]), collapse = "  "), "\n")
    } else {
      for (k in names(tb)) {
        flag <- if (k %in% RESERVED) "   <-- reserved-looking code" else ""
        say(sprintf("      %-10s %6d  (%5.2f%% of observed)%s\n",
                    k, tb[[k]], 100 * tb[[k]] / sum(seen), flag))
      }
    }
  }
  say("\n")
}

say("\n"); line("#")
say("TIME-VARYING VARIABLES\n"); line("#")
for (v in TV) {
  w <- if (v == "sex_birth") 1 else if (v == "sex_current") c(1, 2) else 0:2
  show_var(v, w)
}

say("\n"); line("#")
say("BASELINE-ONLY VARIABLES\n"); line("#")
for (v in TI) show_var(v, 0, invariant = TRUE)

# ---------------------------------------------------------------------------
line("="); say("OUTCOME OBSERVATION PATTERN\n"); line("=")
Y <- c("walk_freq","lsport_freq","msport_freq","ssport_freq")
say("\nhow many of the four frequency outcomes carry a value, per person-wave,\n")
say("counting only person-waves in which the person was observed:\n\n")
for (w in 0:2) {
  d <- long[long$wave == w, ]
  seen <- in_wave[[as.character(w)]][match(d$entity_id, wide$entity_id)]
  d <- d[seen, ]
  k <- rowSums(!is.na(d[, Y]))
  tb <- table(factor(k, levels = 0:4))
  say(sprintf("  wave %d (n = %d):  %s\n", w, nrow(d),
              paste(sprintf("%d of 4 = %d (%.1f%%)", 0:4, as.integer(tb),
                            100 * as.integer(tb) / nrow(d)), collapse = "   ")))
}

say("\nnumber of waves in which a person has all four frequency outcomes present:\n\n")
cnt <- sapply(0:2, function(w) {
  d <- long[long$wave == w, ]
  ok <- rowSums(!is.na(d[, Y])) == 4L
  setNames(ok, d$entity_id)[as.character(wide$entity_id)]
})
nw <- rowSums(cnt, na.rm = TRUE)
tb <- table(factor(nw, levels = 0:3))
for (k in names(tb))
  say(sprintf("  %s wave(s): %6d persons (%.1f%%)\n", k, tb[[k]], 100 * tb[[k]] / nrow(wide)))

say("\nNOTE: 'present' here means the field is not NA. It does NOT yet mean the\n")
say("value is usable -- reserved codes are still counted as present. Deciding\n")
say("which codes become NA is the next step, variable by variable.\n\n")

close(con)
cat(sprintf("\nwritten: %s\n", file.path(OUT, "variable_inventory.txt")))
