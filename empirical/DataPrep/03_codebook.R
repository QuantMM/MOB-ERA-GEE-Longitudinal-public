# ---------------------------------------------------------------------------
# 03_codebook.R
#
# Pull the official CLSA code -> label mapping out of the -dictionary.xlsx
# files and attach it to the observed value counts, so that every code in the
# integrated dataset can be read with its meaning next to it.
#
# Nothing is recoded. This only makes the inventory readable.
#
# Run from the project root:
#     Rscript DataPrep/03_codebook.R
# Writes DataPrep/out/codebook.txt
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(readxl); library(dplyr)})

DATA <- "CLSA Datafiles for analysis"
OUT  <- "DataPrep/out"

DICT <- c(
  w0 = file.path(DATA, "2310007_UofManitoba_SKim_Baseline/2310007_UofManitoba_SKim_Baseline",
                 "2310007_UofManitoba_SKim_Baseline_CoPv7_Qx_CANUE_PA-dictionary.xlsx"),
  w1 = file.path(DATA, "2310007_UofManitoba_SKim_FUP1/2310007_UofManitoba_SKim_FUP1",
                 "2310007_UofManitoba_SKim_FUP1_CoPv5_Qx_CANUE_PA-dictionary.xlsx"),
  w2 = file.path(DATA, "2310007_UofManitoba_SKim_FUP2/2310007_UofManitoba_SKim_FUP2",
                 "2310007_UofManitoba_SKim_FUP2_CoPv2_Qx_PA-dictionary.xlsx")
)

# our name -> source variable name, per wave
MAP <- list(
  walk_freq   = c(w0 = "PA2_WALK_MCQ",     w1 = "PA2_WALK_COF1",     w2 = "PA2_WALK_COF2"),
  walk_hr     = c(w0 = "PA2_WALKHR_MCQ",   w1 = "PA2_WALKHR_COF1",   w2 = "PA2_WALKHR_COF2"),
  lsport_freq = c(w0 = "PA2_LSPRT_MCQ",    w1 = "PA2_LSPRT_COF1",    w2 = "PA2_LSPRT_COF2"),
  lsport_hr   = c(w0 = "PA2_LSPRTHR_MCQ",  w1 = "PA2_LSRTHR_COF1",   w2 = "PA2_LSRTHR_COF2"),
  msport_freq = c(w0 = "PA2_MSPRT_MCQ",    w1 = "PA2_MSPRT_COF1",    w2 = "PA2_MSPRT_COF2"),
  msport_hr   = c(w0 = "PA2_MSPRTHR_MCQ",  w1 = "PA2_MSPRTHR_COF1",  w2 = "PA2_MSPRTHR_COF2"),
  ssport_freq = c(w0 = "PA2_SSPRT_MCQ",    w1 = "PA2_SSPRT_COF1",    w2 = "PA2_SSPRT_COF2"),
  ssport_hr   = c(w0 = "PA2_SSPRTHR_MCQ",  w1 = "PA2_SSPRTHR_COF1",  w2 = "PA2_SSPRTHR_COF2"),
  own         = c(w0 = "OWN_OWN_COM",      w1 = "OWN_OWN_COF1",      w2 = "OWN_OWN_COF2"),
  incneeds    = c(w0 = "WEA_INCNEEDS_MCQ", w1 = "WEA_INCNEEDS_COF1", w2 = "WEA_INCNEEDS_COF2"),
  urban_rural = c(w0 = "SDC_URBAN_RURAL_COM", w1 = "SDC_URBAN_RURAL_COF1",
                  w2 = "SDC_URBAN_RURAL_COF2"),
  age         = c(w0 = "AGE_NMBR_COM",     w1 = "AGE_NMBR_COF1",     w2 = "AGE_NMBR_COF2"),
  sex_current = c(                          w1 = "SDC_CURRSEX_COF1", w2 = "SDC_CURRSEX_COF2"),
  sex_birth   = c(                          w1 = "SDC_BTHSEX_COF1"),
  sex_ask     = c(w0 = "SEX_ASK_COM"),
  yrs_in_canada = c(w0 = "SDC_DRES_COM"),
  cultural_bg = c(w0 = "SDC_DCGT_COM"),
  educ4       = c(w0 = "ED_UDR04_COM"),
  educ_high   = c(w0 = "ED_HIGH_COM"),
  country_birth = c(w0 = "SDC_COB_COM"),
  age_grp     = c(w0 = "AGE_GRP_COM")
)

# --- read the dictionaries --------------------------------------------------
read_dict <- function(path) {
  vars <- read_excel(path, sheet = "Variables", .name_repair = "minimal")
  cats <- read_excel(path, sheet = "Categories", .name_repair = "minimal")
  names(vars) <- tolower(trimws(names(vars)))
  names(cats) <- tolower(trimws(names(cats)))
  list(vars = vars, cats = cats)
}
D <- lapply(DICT, read_dict)

pick <- function(df, want) {
  hit <- names(df)[names(df) %in% want]
  if (length(hit)) df[[hit[1]]] else rep(NA_character_, nrow(df))
}

var_label <- function(w, v) {
  d <- D[[w]]$vars
  nm <- pick(d, c("name", "variable", "variable name"))
  lb <- pick(d, c("label", "description", "variable label"))
  i <- match(v, nm)
  if (is.na(i)) NA_character_ else as.character(lb[i])
}
# NOTE on the Categories sheet layout:
#   columns are  table | variable | name | code | missing | label
#   the CODE VALUE lives in `name` (the `code` column is empty), and `missing`
#   is CLSA's own flag marking which codes are reserved missing values.
code_labels <- function(w, v) {
  d <- D[[w]]$cats
  k <- which(trimws(as.character(d$variable)) == v)
  if (!length(k)) return(NULL)
  lab <- as.character(d$label[k])
  mis <- suppressWarnings(as.integer(d$missing[k]))
  mis[is.na(mis)] <- 0L
  data.frame(code = trimws(as.character(d$name[k])),
             label = lab, missing = mis, stringsAsFactors = FALSE)
}

long <- readRDS(file.path(OUT, "clsa_integrated_long.rds"))
wide <- readRDS(file.path(OUT, "clsa_integrated_wide.rds"))
in_wave <- list(
  "0" = wide$entity_id %in% wide$entity_id,
  "1" = Reduce(`|`, lapply(grep("_w1$", names(wide), value = TRUE), function(v) !is.na(wide[[v]]))),
  "2" = Reduce(`|`, lapply(grep("_w2$", names(wide), value = TRUE), function(v) !is.na(wide[[v]]))))

con <- file(file.path(OUT, "codebook.txt"), open = "wt")
say <- function(...) { t <- paste0(...); cat(t); writeLines(sub("\n$", "", t), con) }
line <- function(ch = "-") say(strrep(ch, 78), "\n")

say("CLSA CODEBOOK + OBSERVED COUNTS\n")
say("official code labels from the -dictionary.xlsx files, next to what the\n")
say("integrated dataset actually contains. Nothing is recoded.\n\n")

for (v in names(MAP)) {
  waves <- names(MAP[[v]])
  line("="); say(sprintf("%s\n", v)); line("=")
  for (w in waves) {
    src <- MAP[[v]][[w]]
    wn  <- as.integer(sub("w", "", w))
    lbl <- var_label(w, src)
    cl  <- code_labels(w, src)

    d <- long[long$wave == wn, ]
    seen <- in_wave[[as.character(wn)]][match(d$entity_id, wide$entity_id)]
    x <- d[[v]]
    tb <- sort(table(x[!is.na(x)]), decreasing = TRUE)

    say(sprintf("\n  wave %d   source variable: %s\n", wn, src))
    if (!is.na(lbl)) say(sprintf("            label: %s\n", lbl))
    say(sprintf("            observed persons %d | item blank %d\n",
                sum(seen), sum(is.na(x) & seen)))
    if (!length(tb)) { say("            (no values)\n"); next }

    labelled <- if (is.null(cl)) character(0) else cl$code
    if (length(tb) > 30) {
      say(sprintf("            %d distinct values; showing every value CLSA gives a\n",
                  length(tb)))
      say("            label to, plus the 8 most frequent\n")
      keep <- union(intersect(names(tb), labelled), names(tb)[1:8])
      tb2 <- tb[names(tb) %in% keep]
    } else tb2 <- tb

    n_mis <- 0L
    for (k in names(tb2)) {
      i <- if (is.null(cl)) NA_integer_ else match(k, cl$code)
      lab <- if (!is.na(i)) cl$label[i] else ""
      flag <- if (!is.na(i) && cl$missing[i] == 1L) " [CLSA: MISSING]" else ""
      if (!is.na(i) && cl$missing[i] == 1L) n_mis <- n_mis + tb2[[k]]
      say(sprintf("      %-8s %7d  %6.2f%%   %s%s\n",
                  k, tb2[[k]], 100 * tb2[[k]] / sum(seen), lab, flag))
    }
    if (length(tb) > length(tb2))
      say(sprintf("      ... %d further values not shown (unlabelled, i.e. genuine data)\n",
                  length(tb) - length(tb2)))
    if (n_mis > 0)
      say(sprintf("      -> flagged missing by CLSA: %d (%.2f%% of observed)\n",
                  n_mis, 100 * n_mis / sum(seen)))
  }
  say("\n")
}
close(con)
cat(sprintf("\nwritten: %s\n", file.path(OUT, "codebook.txt")))
