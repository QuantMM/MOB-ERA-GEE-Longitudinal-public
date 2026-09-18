# ---------------------------------------------------------------------------
# 22_mob_alpha01_sensitivity.R
#
# Sensitivity of the MOB trees to a stricter significance-test stopping rule.
# Every tree in Analysis/16 (D1 full sample, D2 immigrants only) and in
# Analysis/19 Part B (within non-immigrant, recent/mid-term, established) is
# refitted with alpha = .01 instead of .05. Everything else is unchanged:
# same data, candidate sets, reference age, minsize = 100 persons, maxdepth = 4,
# Bonferroni adjustment over candidates, alternating global (gamma, delta)
# estimation.
#
#   PART 1  alpha = .05 trees (saved fits): split variable and adjusted p at
#           every inner node, i.e. which splits a .01 threshold would stop.
#   PART 2  refit at alpha = .01: tree size, splits, inner-node p values,
#           terminal-node estimates.
#
# Root-node tests do not depend on alpha. Counts below 6 are suppressed.
# Run from the project root:
#     Rscript Analysis/22_mob_alpha01_sensitivity.R
# Writes (local only): DataPrep/out/MOB_ALPHA01_SENSITIVITY.md, _mob_alpha01.rds
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("FINAL_MOB-ERA-GEE_Code/Empirical Application/FINAL_HELPERS_MOB_ERA_GEE.R")
source("Analysis/ERA_GEE_fixed.R")
source("Analysis/nuisance_adjusted_slope.R")
source("Analysis/palm_era_gee.R")

OUT <- "DataPrep/out"
Y   <- c("y_walk", "y_lsport", "y_msport", "y_ssport")
YL  <- c("walk", "light", "mod", "stren")
Q   <- length(Y); FAM <- rep("gaussian", Q)
MINSIZE <- 100L; MAXDEPTH <- 4L; ALPHA <- 0.01

sup <- function(n) if (n > 0 && n < 6) "<6" else format(n, big.mark = ",", trim = TRUE)

# --- data, exactly as in Analysis/16 and Analysis/19 ------------------------
d   <- readRDS(file.path(OUT, "analytic_long_balanced.rds"))
raw <- readRDS(file.path(OUT, "immigration_items_raw.rds"))
i   <- match(d$entity_id, raw$entity_id)
fl  <- raw$SDC_FIMM_COM[i]; yrs <- suppressWarnings(as.numeric(raw$SDC_DRES_COM[i]))
d$z_imm_flag <- factor(ifelse(fl == "1", "Immigrant", ifelse(fl == "2", "Non-immigrant", NA)),
                       levels = c("Non-immigrant", "Immigrant"))
d$z_duration <- factor(ifelse(is.na(yrs) | yrs >= 996, NA,
                       ifelse(yrs < 20, "Recent/mid-term", "Long-term")),
                       levels = c("Long-term", "Recent/mid-term"))
d$grp <- factor(ifelse(fl == "2", "Non-immigrant",
                ifelse(fl == "1" & !is.na(yrs) & yrs < 996,
                       ifelse(yrs < 20, "Recent/mid-term", "Established"), NA)),
                levels = c("Non-immigrant", "Recent/mid-term", "Established"))
AGE_C <- mean(d$z_age[d$wave == 0])

Z_COMMON <- c("z_sex", "z_racialized", "z_educ4", "z_incneeds", "z_own")
SETS <- list(
  `D1`              = list(data = d,                                                    Z = c("z_imm_flag", Z_COMMON)),
  `D2`              = list(data = d[!is.na(d$z_imm_flag) & d$z_imm_flag == "Immigrant", ], Z = c("z_duration", Z_COMMON)),
  `Non-immigrant`   = list(data = d[!is.na(d$grp) & d$grp == "Non-immigrant", ],       Z = Z_COMMON),
  `Recent/mid-term` = list(data = d[!is.na(d$grp) & d$grp == "Recent/mid-term", ],     Z = Z_COMMON),
  `Established`     = list(data = d[!is.na(d$grp) & d$grp == "Established", ],         Z = Z_COMMON))
DG <- c("joint", "slope", "intercept")

# --- helpers ----------------------------------------------------------------
inner_nodes <- function(tr) {
  ids <- setdiff(nodeids(tr), nodeids(tr, terminal = TRUE))
  if (!length(ids)) return("no split")
  paste(sapply(ids, function(k) {
    nd <- nodeapply(tr, k, identity)[[1]]; ti <- info_node(nd); sp <- split_node(nd)
    vn <- names(tr$data)[sp$varid]; zz <- tr$data[[sp$varid]]
    lab <- if (is.factor(zz)) paste0("{", paste(levels(zz)[which(sp$index == 1L)], collapse = ","), "}|{",
                                     paste(levels(zz)[which(sp$index == 2L)], collapse = ","), "}") else ""
    sprintf("node %d [%s persons] %s %s p=%.3g", k, sup(ti$nobs %/% 3), vn, lab, ti$test["p.value", vn])
  }), collapse = "\n      ")
}
terminal_report <- function(fit, dat) {
  node <- fit$node; off <- fit$offset; out <- character()
  for (nd in sort(unique(node[!is.na(node)]))) {
    ix <- which(!is.na(node) & node == nd); s <- dat[ix, ]; np <- length(unique(s$id))
    if (np < 6) { out <- c(out, sprintf("node %d: [suppressed]", nd)); next }
    m <- ERA_GEE_fixed(as.matrix(s[, Y]) - off[ix, , drop = FALSE], s$X_time, FAM, "independence", y_scale = "raw")
    out <- c(out, sprintf("node %-2d n=%7s  %%F=%3.0f %%rac=%3.0f  b1 %s | b2 %s", nd, sup(np),
                          100 * mean(s$z_sex == "Female"), 100 * mean(s$z_racialized == "Racialized", na.rm = TRUE),
                          paste(sprintf("%.3f", m$B["intercept", ]), collapse = " "),
                          paste(sprintf("%+.3f", m$B["slope", ]), collapse = " ")))
  }
  paste(out, collapse = "\n      ")
}

con <- file(file.path(OUT, "MOB_ALPHA01_SENSITIVITY.md"), open = "wt")
w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con); cat(s, "\n") }
w("# MOB sensitivity: alpha = .01 versus alpha = .05"); w("")
w("**Generated:** ", format(Sys.time(), "%d %B %Y %H:%M"), " - `Analysis/22_mob_alpha01_sensitivity.R`"); w("")
w("Outcome order in node estimates: ", paste(YL, collapse = ", "), ". b1 = fitted Wave-1 level, b2 = change per wave, at age ", sprintf("%.2f", AGE_C), "."); w("")

# ===========================================================================
w("## PART 1: saved alpha = .05 trees, inner-node split p values"); w("")
m05 <- readRDS(file.path(OUT, "_mob_final_model.rds"))
t05 <- readRDS(file.path(OUT, "_three_group.rds"))$mob
old <- list()
for (dg in DG) { old[[paste("D1", dg)]] <- m05[[paste("D1", dg)]]; old[[paste("D2", dg)]] <- m05[[paste("D2", dg)]] }
for (g in c("Non-immigrant", "Recent/mid-term", "Established")) for (dg in DG) old[[paste(g, dg)]] <- t05[[paste(g, dg)]]
for (k in names(old)) {
  tr <- old[[k]]$tree
  w("- **", k, "**: ", width(tr), " terminal nodes"); w("      ", inner_nodes(tr))
}

# ===========================================================================
w(""); w("## PART 2: refit at alpha = .01"); w("")
ctrl <- function(prm) mob_control(verbose = FALSE, maxdepth = MAXDEPTH, minsize = MINSIZE * 3L,
                                  alpha = ALPHA, parm = prm)
res <- list()
for (sn in names(SETS)) {
  dd <- droplevels(SETS[[sn]]$data); dd$rowid <- seq_len(nrow(dd))
  for (dg in DG) {
    prm <- switch(dg, joint = 1:(2 * Q), intercept = 1:Q, slope = NULL)
    t0 <- Sys.time()
    f <- palm_era_gee_mob(dd, SETS[[sn]]$Z, dg, Y, "X_time", FAM, age_var = "z_age",
                          age_center = AGE_C, control = ctrl(prm))
    k <- paste(sn, dg); res[[k]] <- f
    w("- **", k, "**: ", width(f$tree), " terminal nodes (alpha .05: ", width(old[[k]]$tree), ")  [",
      sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ", ",
      if (isTRUE(f$converged)) "converged" else "NOT CONVERGED", "]")
    w("      ", inner_nodes(f$tree))
    if (width(f$tree) > 1L) w("      ", terminal_report(f, dd))
  }
}
saveRDS(lapply(res, function(f) f[c("node", "gamma", "converged", "iter")]), file.path(OUT, "_mob_alpha01.rds"))
close(con)
cat("\nwritten: DataPrep/out/MOB_ALPHA01_SENSITIVITY.md\n")
