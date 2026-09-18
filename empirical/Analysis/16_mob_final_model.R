# ---------------------------------------------------------------------------
# 16_mob_final_model.R
#
# MOB on top of the FINAL foundational ERA-GEE model.
#
# MEAN MODEL (identical to the foundational model, Analysis/15):
#     g(mu_itq) = beta_1q + beta_2q t + gamma_q A_c + delta_q (t x A_c)
#
#   gamma_q, delta_q  GLOBAL: estimated once, common to every terminal node
#   beta_1q, beta_2q  NODE-SPECIFIC: these are the parameters MOB tests
#
# Age is therefore an adjustment carried in the mean model, NOT a partitioning
# variable. The question is what trajectory heterogeneity remains once level
# and change differences attributable to baseline-visit age are removed.
# Estimation uses the alternating scheme in Analysis/palm_era_gee.R.
#
# TWO DESIGNS
#   D1  full sample; immigrant status from the CLSA Immigration Flag
#       (SDC_FIMM_COM: 1 = immigrant, 2 = non-immigrant, 9 = not answered)
#   D2  immigrants only; settlement-duration category
#       (recent/mid-term < 20 years vs long-term >= 20 years)
#
#   Both use the same remaining partitioning variables: sex, racialization,
#   education, income adequacy, and dwelling ownership. Area of residence is
#   not included. Racialization enters as the binary contrast; the
#   disaggregated population-group variable is never entered alongside it.
#
# Trees are grown permissively and ARE NOT PRUNED: minsize = 100 persons,
# maxdepth = 4, alpha = .05 with Bonferroni adjustment over the partitioning
# variables. Counts below 6 are suppressed.
#
# Run from the project root:
#     Rscript Analysis/16_mob_final_model.R
# Writes (local only): DataPrep/out/MOB_FINAL_MODEL.md
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("FINAL_MOB-ERA-GEE_Code/Empirical Application/FINAL_HELPERS_MOB_ERA_GEE.R")
source("Analysis/ERA_GEE_fixed.R")
source("Analysis/nuisance_adjusted_slope.R")
source("Analysis/palm_era_gee.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c("walk", "lsport", "msport", "ssport")
Q   <- length(Y); FAM <- rep("gaussian", Q)
MINSIZE_PERSONS <- 100L; MAXDEPTH <- 4L; ALPHA <- 0.05

sup  <- function(n) if (n > 0 && n < 6) "<6" else format(n, big.mark = ",", trim = TRUE)
line <- function(ch = "-") cat(strrep(ch, 78), "\n")

# ---------------------------------------------------------------------------
# 1. Data: balanced sample + the two immigration specifications
# ---------------------------------------------------------------------------
d   <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
raw <- readRDS(file.path(OUT, "immigration_items_raw.rds"))

fimm <- raw$SDC_FIMM_COM[match(d$entity_id, raw$entity_id)]
d$z_imm_flag <- factor(ifelse(fimm == "1", "Immigrant",
                       ifelse(fimm == "2", "Non-immigrant", NA)),
                       levels = c("Non-immigrant", "Immigrant"))

yrs <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[match(d$entity_id, raw$entity_id)]))
d$z_duration <- factor(ifelse(is.na(yrs) | yrs >= 996, NA,
                       ifelse(yrs < 20, "Recent/mid-term", "Long-term")),
                       levels = c("Long-term", "Recent/mid-term"))

AGE_C <- mean(d$z_age[d$wave == 0])     # same reference age as Analysis/15

Z_COMMON <- c("z_sex", "z_racialized", "z_educ4", "z_incneeds", "z_own")
DESIGNS <- list(
  D1 = list(lab  = "Full sample; immigrant status from the CLSA Immigration Flag",
            sub  = FALSE,
            Z    = c("z_imm_flag", Z_COMMON)),
  D2 = list(lab  = "Immigrants only; recent/mid-term vs long-term settlement",
            sub  = TRUE,
            Z    = c("z_duration", Z_COMMON))
)

ctrl <- function(prm) mob_control(verbose = FALSE, maxdepth = MAXDEPTH,
                                  minsize = MINSIZE_PERSONS * 3L,
                                  alpha = ALPHA, parm = prm)

# ---------------------------------------------------------------------------
# 2. Reporting helpers
# ---------------------------------------------------------------------------
split_desc <- function(tr) {
  inner <- setdiff(nodeids(tr), nodeids(tr, terminal = TRUE))
  if (!length(inner)) return("(no split)")
  paste(sapply(inner, function(i) {
    nd <- nodeapply(tr, i, identity)[[1]]; sp <- split_node(nd)
    vn <- names(tr$data)[sp$varid]; zz <- tr$data[[sp$varid]]
    if (is.factor(zz)) {
      k <- sp$index
      sprintf("node %d: %s  {%s} | {%s}", i, vn,
              paste(levels(zz)[which(k == 1L)], collapse = ","),
              paste(levels(zz)[which(k == 2L)], collapse = ","))
    } else sprintf("node %d: %s <= %.4g", i, vn, sp$breaks[1])
  }), collapse = "  ||  ")
}

node_report <- function(fit, dat, zvars) {
  tr <- fit$tree; node <- fit$node; off <- fit$offset
  cat(sprintf("\n  terminal nodes: per-wave change b2 (top) and fitted Wave-1 level b1 (bottom),\n"))
  cat(sprintf("  both at baseline-visit age %.2f\n", AGE_C))
  cat(sprintf("  %-5s %8s %6s %5s %5s %5s  %s\n", "node", "n", "age", "%F", "%imm", "%rac",
              paste(sprintf("%8s", YL), collapse = " ")))
  for (nd in sort(unique(node[!is.na(node)]))) {
    ix <- which(!is.na(node) & node == nd)
    dd <- dat[ix, ]; np <- length(unique(dd$id))
    if (np < 6) { cat(sprintf("  %-5d %8s  [suppressed, n < 6]\n", nd, "<6")); next }
    yy <- as.matrix(dd[, Y, drop = FALSE]) - off[ix, , drop = FALSE]
    m  <- ERA_GEE_fixed(yy, dd$X_time, FAM, "independence", y_scale = "raw")
    pim <- if ("z_imm_flag" %in% names(dd)) 100 * mean(dd$z_imm_flag == "Immigrant", na.rm = TRUE) else NA
    cat(sprintf("  %-5d %8s %6.1f %5.0f %5.0f %5.0f  %s\n", nd, format(np, big.mark = ","),
                mean(dd$z_age), 100 * mean(dd$z_sex == "Female", na.rm = TRUE), pim,
                100 * mean(dd$z_racialized == "Racialized", na.rm = TRUE),
                paste(sprintf("%+8.4f", m$B["slope", ]), collapse = " ")))
    cat(sprintf("  %-5s %8s %6s %5s %5s %5s  %s\n", "", "", "", "", "", "",
                paste(sprintf("%8.4f", m$B["intercept", ]), collapse = " ")))
  }
  # composition of the split variables actually used
  used <- unique(unlist(lapply(setdiff(nodeids(tr), nodeids(tr, terminal = TRUE)),
    function(i) names(tr$data)[split_node(nodeapply(tr, i, identity)[[1]])$varid])))
  if (length(used)) cat("\n  split variables used:", paste(used, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# 3. Fit
# ---------------------------------------------------------------------------
res <- list()
for (nm in names(DESIGNS)) {
  s  <- DESIGNS[[nm]]
  dd <- d
  if (s$sub) dd <- dd[!is.na(dd$z_imm_flag) & dd$z_imm_flag == "Immigrant", ]
  dd <- droplevels(dd); dd$rowid <- seq_len(nrow(dd))

  line("#"); cat(sprintf("DESIGN %s -- %s\n", nm, s$lab)); line("#")
  cat(sprintf("  %s persons, %s person-waves | Z = %s\n",
              format(length(unique(dd$id)), big.mark = ","),
              format(nrow(dd), big.mark = ","), paste(s$Z, collapse = ", ")))
  nmiss <- sapply(s$Z, function(v) sum(is.na(dd[[v]][!duplicated(dd$id)])))
  cat("  persons with a missing value on each Z: ",
      paste(sprintf("%s=%s", names(nmiss), sapply(nmiss, sup)), collapse = ", "), "\n")
  cat(sprintf("  mob() drops person-waves missing on ANY Z: %s of %s persons\n",
              sup(sum(!complete.cases(dd[!duplicated(dd$id), s$Z]))),
              format(length(unique(dd$id)), big.mark = ",")))
  for (v in s$Z) {
    tb <- table(dd[[v]][!duplicated(dd$id)])
    cat(sprintf("    %-14s %s\n", v, paste(sprintf("%s=%s", names(tb), sapply(tb, sup)), collapse = "  ")))
  }

  for (dg in c("joint", "slope", "intercept")) {
    prm <- switch(dg, joint = 1:(2 * Q), intercept = 1:Q, slope = NULL)
    t0  <- Sys.time()
    fit <- palm_era_gee_mob(dd, s$Z, dg, Y, "X_time", FAM,
                            age_var = "z_age", age_center = AGE_C, control = ctrl(prm))
    el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    tr  <- fit$tree
    ti  <- nodeapply(tr, 1, function(n) info_node(n))[[1]]
    o   <- order(ti$test["p.value", ])

    cat("\n"); line("-")
    cat(sprintf("%s | %s diagnostic   (%.1f min, %s in %d iterations)\n", nm, dg, el,
                if (fit$converged) "converged" else "NOT CONVERGED", fit$iter))
    line("-")
    cat("  root parameter instability (sorted by p):\n")
    for (j in o) cat(sprintf("    %-14s statistic = %8.2f   p = %.3g\n",
                             colnames(ti$test)[j], ti$test["statistic", j], ti$test["p.value", j]))
    cat(sprintf("\n  terminal nodes: %d\n  splits: %s\n", width(tr), split_desc(tr)))
    cat("\n  global age effects (per year of baseline-visit age):\n")
    cat(sprintf("    %-16s %s\n", "", paste(sprintf("%9s", YL), collapse = " ")))
    cat(sprintf("    %-16s %s\n", "gamma (level)",
                paste(sprintf("%+9.5f", fit$gamma["age", ]), collapse = " ")))
    cat(sprintf("    %-16s %s\n", "delta (per wave)",
                paste(sprintf("%+9.5f", fit$gamma["age_time", ]), collapse = " ")))
    if (width(tr) > 1L) node_report(fit, dd, s$Z)
    res[[paste(nm, dg)]] <- fit
  }
}

saveRDS(res, file.path(OUT, "_mob_final_model.rds"))
cat("\n"); line("="); cat("saved: DataPrep/out/_mob_final_model.rds\n"); line("=")
