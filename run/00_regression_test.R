# ---------------------------------------------------------------------------
# 00c_regression_test.R
#
# REGRESSION TEST for the reference case.
#
# Run this after ANY change to ERA_GEE_fixed.R, synthetic_data.R, or the
# patches. It asserts that the frozen reference example still recovers the
# intended answer:
#
#     IMM is the first split, cutting {recent} from {non, established},
#     with terminal-node intercepts and slopes matching the known truth.
#
# Exits non-zero on failure.
#     Rscript run/00_regression_test.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(partykit); library(MASS)})
source("R/era_gee_fixed.R")
source("R/synthetic_data.R")

ok <- TRUE
check <- function(label, cond, detail = "") {
  cond <- isTRUE(cond)
  cat(sprintf("  [%s] %s%s\n", if (cond) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  --  ", detail) else ""))
  if (!cond) ok <<- FALSE
  invisible(cond)
}

ref <- make_reference_data()
dat <- ref$dat; YN <- ref$YN; ZN <- ref$ZN; Q <- ref$Q
fam <- rep("gaussian", Q)

cat("REGRESSION TEST -- reference case\n")
cat(sprintf("  data: %d persons x %d waves = %d rows, Q = %d\n\n",
            ref$n_person, ref$n_wave, nrow(dat), Q))

# --- 1. root model shape --------------------------------------------------
cat("1. root model\n")
root <- ERA_GEE_fixed(dat[, YN], dat$X_time, fam, "independence", y_scale = "raw")
check("B is 2 x Q",            identical(dim(root$B), c(2L, as.integer(Q))))
check("rownames of B",         identical(rownames(root$B), c("intercept", "slope")))
check("estfun is n_row x 2Q",  identical(dim(root$estfun), c(nrow(dat), 2L * as.integer(Q))))
check("scores vanish at optimum", max(abs(colSums(root$estfun))) < 1e-6,
      sprintf("max |colSum| = %.2e", max(abs(colSums(root$estfun)))))
check("f_const == 1 exactly",  max(abs(root$F[, "f_const"] - 1)) == 0)
check("f_time  == t exactly",  max(abs(root$F[, "f_time"] - dat$X_time)) == 0)

# --- 2. tree structure ----------------------------------------------------
cat("\n2. MOB tree (method A: wave-level estfun + cluster = id)\n")
fitf <- era_gee_mob_fit(dat, YN, "X_time", fam, "independence", y_scale = "raw")
fml  <- as.formula(paste0("cbind(", paste(YN, collapse = ","), ") ~ rowid + X_time | ",
                          paste(ZN, collapse = " + ")))
tr <- mob(fml, data = dat, fit = fitf, cluster = dat$id,
          control = mob_control(verbose = FALSE, maxdepth = 3, minsize = 450))

check("exactly 2 terminal nodes", width(tr) == 2L, sprintf("got %d", width(tr)))

# NOTE: split_node()$varid indexes the mob MODEL FRAME (response, regressors,
# then partitioning variables), not the original data frame.
sv <- names(tr$data)[nodeapply(tr, 1, function(n) split_node(n)$varid)[[1]]]
check("first split variable is IMM", identical(sv, "IMM"), sprintf("got '%s'", sv))

memb   <- predict(tr, type = "node")
leaves <- nodeids(tr, terminal = TRUE)
grp <- lapply(leaves, function(nd) {
  d <- dat[memb == nd, ]; sort(unique(as.character(d$IMM[!duplicated(d$id)])))
})
check("one node is exactly {recent}",
      any(sapply(grp, function(g) identical(g, "recent"))))
check("other node is exactly {established, non}",
      any(sapply(grp, function(g) identical(g, c("established", "non")))))

# --- 3. node parameters vs truth ------------------------------------------
cat("\n3. terminal-node parameters vs truth (tol: intercept 0.15, slope 0.05)\n")
for (k in seq_along(leaves)) {
  nd <- leaves[k]; d <- dat[memb == nd, ]
  m  <- ERA_GEE_fixed(d[, YN], d$X_time, fam, "independence", y_scale = "raw")
  is_rec  <- identical(grp[[k]], "recent")
  tb1 <- if (is_rec) ref$truth$b1_recent else ref$truth$b1
  tb2 <- if (is_rec) ref$truth$b2_recent else ref$truth$b2_others
  di  <- max(abs(m$B["intercept", ] - tb1))
  ds  <- max(abs(m$B["slope", ]     - tb2))
  check(sprintf("node %d {%s}: intercepts", nd, paste(grp[[k]], collapse = "+")),
        di < 0.15, sprintf("max abs dev = %.3f", di))
  check(sprintf("node %d {%s}: slopes",     nd, paste(grp[[k]], collapse = "+")),
        ds < 0.05, sprintf("max abs dev = %.3f", ds))
}

cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASSED" else "*** FAILURES ABOVE ***"))
if (!ok) quit(status = 1L)
