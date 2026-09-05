# ---------------------------------------------------------------------------
# 06_adopted_diagnostics_table.R
#
# The 4 x 3 table under the ADOPTED diagnostic definitions, plus Option A as a
# comparison column. This is the artefact the repeated simulation is built on.
#
# ADOPTED (all three on raw t = 0,1,2, all with cluster = id):
#   joint      full 2Q model, test all 2Q parameters      parm = 1:(2Q)
#   slope      full 2Q model, test the Q nuisance-        no parm; the adjusted
#              adjusted slope scores U_S.I               Q-dim process is passed
#   intercept  full 2Q model, test the Q intercepts       parm = 1:Q
#
# COMPARISON (not adopted, kept as a cross-check):
#   slopeA     full 2Q model on CENTRED t, naive parm = (Q+1):(2Q)
#
# EXPECTED for the adopted slope column (fixed in advance):
#   slope-only -> IMM,  intercept-only -> (none),
#   intercept+slope -> IMM,  null -> (none)
#
# Run from the project root:
#     Rscript run/04_diagnostics_table.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("R/era_gee_fixed.R")
source("R/nuisance_adjusted_slope.R")
source("R/synthetic_data.R")

line <- function(ch = "-") cat(strrep(ch, 78), "\n")

B1     <- REF_TRUTH$b1
B2_OTH <- REF_TRUTH$b2_others
B2_REC <- REF_TRUTH$b2_recent
B1_REC <- B1 - c(0.40, 0.30, 0.25, 0.35)

CONDITIONS <- list(
  `slope-only`      = list(b1_recent = B1,     b2_recent = B2_REC),
  `intercept-only`  = list(b1_recent = B1_REC, b2_recent = B2_OTH),
  `intercept+slope` = list(b1_recent = B1_REC, b2_recent = B2_REC),
  `null`            = list(b1_recent = B1,     b2_recent = B2_OTH)
)

EXPECTED <- rbind(`slope-only`      = c("IMM",    "IMM",    "(none)"),
                  `intercept-only`  = c("IMM",    "(none)", "IMM"),
                  `intercept+slope` = c("IMM",    "IMM",    "IMM"),
                  `null`            = c("(none)", "(none)", "(none)"))
colnames(EXPECTED) <- c("joint", "slope", "intercept")

MINSIZE_ROW <- 450L; MAXDEPTH <- 3L

run_tree <- function(d, YN, ZN, Q, fitf, parm = NULL) {
  fml <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                           paste(ZN, collapse = " + ")))
  tr <- mob(fml, data = d, fit = fitf, cluster = d$id,
            control = mob_control(verbose = FALSE, maxdepth = MAXDEPTH,
                                  minsize = MINSIZE_ROW, parm = parm))
  ti <- nodeapply(tr, 1, function(n) info_node(n))[[1]]
  sv <- if (width(tr) == 1L) "(none)"
        else names(tr$data)[nodeapply(tr, 1, function(n) split_node(n)$varid)[[1]]]
  gp <- if (width(tr) == 1L) "-" else {
    mb <- predict(tr, type = "node")
    paste(sapply(nodeids(tr, terminal = TRUE), function(nd) {
      dd <- d[mb == nd, ]
      paste(sort(unique(as.character(dd$IMM[!duplicated(dd$id)]))), collapse = "+")
    }), collapse = " | ")
  }
  list(split = sv, groups = gp, nnodes = width(tr),
       stat = ti$test["statistic", "IMM"], p = ti$test["p.value", "IMM"])
}

rows <- list()

for (cn in names(CONDITIONS)) {
  cc  <- CONDITIONS[[cn]]
  ref <- make_reference_data(b1_recent = cc$b1_recent, b2_recent = cc$b2_recent)
  dat <- ref$dat; YN <- ref$YN; ZN <- ref$ZN; Q <- ref$Q
  fam <- rep("gaussian", Q)
  datc <- dat; datc$X_time <- datc$X_time - mean(sort(unique(datc$X_time)))

  line("="); cat(sprintf("CONDITION: %s\n", toupper(cn))); line("=")
  print(round(rbind(`intercept: others` = B1, `intercept: recent` = cc$b1_recent,
                    `slope: others` = B2_OTH, `slope: recent` = cc$b2_recent), 3))
  cat("\n")

  specs <- list(
    joint     = list(d = dat,  f = era_gee_mob_fit(dat, YN, "X_time", fam,
                                                   "independence", y_scale = "raw"),
                     parm = 1:(2 * Q), adopted = TRUE),
    slope     = list(d = dat,  f = era_gee_mob_fit_slopeadj(dat, YN, "X_time", fam),
                     parm = NULL, adopted = TRUE),
    intercept = list(d = dat,  f = era_gee_mob_fit(dat, YN, "X_time", fam,
                                                   "independence", y_scale = "raw"),
                     parm = 1:Q, adopted = TRUE),
    slopeA    = list(d = datc, f = era_gee_mob_fit(datc, YN, "X_time", fam,
                                                   "independence", y_scale = "raw"),
                     parm = (Q + 1):(2 * Q), adopted = FALSE)
  )

  for (sn in names(specs)) {
    s  <- specs[[sn]]
    t0 <- Sys.time()
    r  <- run_tree(s$d, YN, ZN, Q, s$f, s$parm)
    ev <- if (s$adopted) EXPECTED[cn, sn] else EXPECTED[cn, "slope"]
    ok <- identical(r$split, ev)
    cat(sprintf("  %-10s%-9s -> %-8s nodes=%d  IMM stat=%9.2f  p=%-11.3g exp %-8s %s (%.0fs)\n",
                sn, if (s$adopted) "[adopted]" else "[compare]",
                r$split, r$nnodes, r$stat, r$p, ev,
                if (ok) "OK" else "<<< MISMATCH",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (r$nnodes > 1L) cat(sprintf("                       terminal groups: %s\n", r$groups))
    rows[[length(rows) + 1L]] <- data.frame(
      condition = cn, diagnostic = sn, adopted = s$adopted, split = r$split,
      expected = ev, match = ok, n_nodes = r$nnodes,
      IMM_stat = round(r$stat, 2), IMM_p = signif(r$p, 3))
  }
  cat("\n")
}

tab <- do.call(rbind, rows)
line("="); cat("ADOPTED 4 x 3 TABLE\n"); line("=")
ad  <- subset(tab, adopted)
got <- with(ad, tapply(split, list(condition, diagnostic), identity))[
         rownames(EXPECTED), colnames(EXPECTED)]
cat("\nOBSERVED (split variable):\n"); print(got)
cat("\nEXPECTED:\n"); print(EXPECTED)
cat(sprintf("\nMATCHES IN ALL 12 ADOPTED CELLS: %s\n", identical(unname(got), unname(EXPECTED))))

cat("\nIMM statistic:\n")
print(with(ad, tapply(IMM_stat, list(condition, diagnostic), identity))[
        rownames(EXPECTED), colnames(EXPECTED)])

line("="); cat("SLOPE DIAGNOSTIC: adopted (Option B) vs comparison (Option A)\n"); line("=")
cmp <- merge(subset(tab, diagnostic == "slope",  c(condition, split, IMM_stat)),
             subset(tab, diagnostic == "slopeA", c(condition, split, IMM_stat)),
             by = "condition", suffixes = c("_B", "_A"))
cmp <- cmp[match(rownames(EXPECTED), cmp$condition), ]
cmp$pct_diff <- round(100 * abs(cmp$IMM_stat_A - cmp$IMM_stat_B) / cmp$IMM_stat_B, 1)
cmp$same_conclusion <- cmp$split_B == cmp$split_A
print(cmp, row.names = FALSE)

saveRDS(tab, "results/_adopted_diagnostics_table.rds")
