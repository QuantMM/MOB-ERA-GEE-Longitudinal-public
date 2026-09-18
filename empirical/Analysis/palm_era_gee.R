# ---------------------------------------------------------------------------
# palm_era_gee.R
#
# ERA-GEE with a GLOBAL age adjustment and a node-specific intercept/slope.
#
#     g(mu_itq) = beta_1q + beta_2q t + gamma_1q Age_i + gamma_2q (Age_i x t)
#
#     f_1 = 1        f_2 = t        f_3 = Age_i        f_4 = Age_i x t
#     all four component weights fixed at 1, exactly as in the base model.
#
# WHAT IS GLOBAL AND WHAT IS LOCAL
#   gamma_1q, gamma_2q   estimated ONCE, common to every terminal node
#   beta_1q,  beta_2q    node-specific; these are the parameters whose
#                        instability MOB tests
#
#   The substantive question this answers is: once differences in level and in
#   rate of change that are attributable to BASELINE AGE have been removed
#   globally, is there any immigration-related trajectory heterogeneity left?
#
# WHY AN ALTERNATING SCHEME IS NEEDED
#   mob() refits its fit-function inside every node, so a parameter that must
#   be COMMON across nodes cannot be estimated by mob() itself. The structure
#   is the partially-additive one that partykit::palmtree implements for
#   linear models, and the algorithm here is the same:
#
#     0. gamma <- 0
#     1. offset  o_itq = gamma_1q Age_i + gamma_2q Age_i t
#     2. grow the MOB tree on (y - o) with the ordinary TWO-component
#        ERA-GEE fitter -- unchanged, including the nuisance-adjusted slope
#        diagnostic
#     3. holding the terminal-node assignment fixed, refit {beta_m, gamma}
#        jointly and take the new gamma
#     4. repeat 1-3 until the node assignment and gamma both stop changing
#
#   Steps 2 and 3 are alternating minimisation of one least-squares objective
#   (Gaussian family, working independence => the ERA-GEE update IS OLS), so
#   each step is exact; only the tree re-growth breaks strict monotonicity,
#   which is why convergence is declared on node assignment as well as gamma.
#
# WHY AGE IS CENTRED
#   With raw age, beta_1q is the level at age 0 and beta_2q the slope at age 0
#   -- extrapolations roughly 45 years outside the data, and nearly collinear
#   with gamma. Centring at the analysis-sample mean baseline age makes
#   beta_1q the level and beta_2q the change at the average age, which is the
#   quantity the instability test should be about.
#
#   Centring does NOT change the test. Under Age -> Age - a the design becomes
#   F A with A = [[1,0,-a,0],[0,1,0,-a],[0,0,1,0],[0,0,0,1]], so the scores
#   transform as A' U, giving U_beta1,new = U_beta1 and U_beta2,new = U_beta2:
#   the beta block is invariant elementwise. Only the gamma-block scores and
#   the reported beta VALUES move. Same argument as for the time origin.
#
# LIMITATION, STATED EXPLICITLY
#   The node-level scores treat gamma as KNOWN. The instability test therefore
#   conditions on gamma-hat and does not propagate its uncertainty. This is the
#   standard partially-additive-tree caveat. It matters most at the root, where
#   the node and the global sample coincide; `pernode_gamma()` below is
#   provided to check the assumption that one gamma fits every node.
# ---------------------------------------------------------------------------

# --- age design: the two additional fixed components -----------------------
make_age_F <- function(age_c, time)
  cbind(age = as.numeric(age_c), age_time = as.numeric(age_c) * as.numeric(time))


# --- step 3: refit {beta_m, gamma} jointly, return gamma (2 x Q) -----------
# Gaussian + working independence, so the ERA-GEE update is exactly OLS of
# each outcome on [node dummies, node dummies x t, Age_c, Age_c x t].
update_gamma <- function(y, time, Fa, node) {
  # built by hand rather than with model.matrix(): a tree that has not split
  # leaves a single-level factor, which model.matrix() refuses to contrast.
  nf <- factor(node)
  D  <- vapply(levels(nf), function(l) as.numeric(nf == l), numeric(length(nf)))
  D  <- matrix(D, nrow = length(nf))
  X <- cbind(D, D * as.numeric(time), Fa)
  colnames(X) <- c(paste0("n", seq_len(ncol(D))),
                   paste0("n", seq_len(ncol(D)), "_t"), colnames(Fa))
  cf <- qr.coef(qr(X), as.matrix(y))
  g  <- cf[colnames(Fa), , drop = FALSE]
  g[is.na(g)] <- 0                      # aliased column -> no contribution
  g
}


# --- diagnostic: does ONE gamma fit every node? ----------------------------
# Unconstrained per-node age effects, fitted after convergence. If these are
# close to the global gamma the parallel-age assumption is doing no damage; if
# they diverge, the global adjustment is misspecified for this partition.
pernode_gamma <- function(y, time, Fa, node) {
  y <- as.matrix(y)
  out <- lapply(sort(unique(node)), function(m) {
    ix <- node == m
    X  <- cbind(`(Intercept)` = 1, t = as.numeric(time)[ix], Fa[ix, , drop = FALSE])
    cf <- qr.coef(qr(X), y[ix, , drop = FALSE])
    cf <- cf[colnames(Fa), , drop = FALSE]
    cf[is.na(cf)] <- 0
    cf
  })
  names(out) <- as.character(sort(unique(node)))
  out
}


# --- mob fit factories, identical to the base ones but on (y - offset) -----
# `off` is an nrow(dat.full) x Q matrix of the current global age contribution.
era_gee_mob_fit_off <- function(dat.full, ynames, time_var, family, off,
                                corstr = "independence", rowid_var = "rowid") {
  force(off)
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) {
    node_rows <- as.integer(round(as.matrix(x)[, rowid_var]))
    d  <- dat.full[node_rows, , drop = FALSE]
    yy <- as.matrix(d[, ynames, drop = FALSE]) - off[node_rows, , drop = FALSE]
    m  <- ERA_GEE_fixed(y = yy, time = d[[time_var]], family = family,
                        corstr = corstr, y_scale = "raw")
    rval <- list(coefficients = setNames(m$coef, colnames(m$estfun)),
                 objfun = m$objfun)
    if (estfun) rval$estfun <- m$estfun
    if (object) rval$object <- m
    rval
  }
}

era_gee_mob_fit_slopeadj_off <- function(dat.full, ynames, time_var, family, off,
                                         corstr = "independence",
                                         rowid_var = "rowid", id_var = "id") {
  force(off); Q <- length(ynames)
  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) {
    node_rows <- as.integer(round(as.matrix(x)[, rowid_var]))
    d  <- dat.full[node_rows, , drop = FALSE]
    yy <- as.matrix(d[, ynames, drop = FALSE]) - off[node_rows, , drop = FALSE]
    m  <- ERA_GEE_fixed(y = yy, time = d[[time_var]], family = family,
                        corstr = corstr, y_scale = "raw")
    adj <- nuisance_adjust_slope(m$estfun, d[[id_var]], Q)
    rval <- list(coefficients = setNames(m$B["slope", ], colnames(adj)),
                 objfun = m$objfun)
    if (estfun) rval$estfun <- adj
    if (object) rval$object <- m
    rval
  }
}


# --- the alternating fit ---------------------------------------------------
# Returns the converged tree, the global gamma, and the per-node gammas.
palm_era_gee_mob <- function(d, zvars, diagnostic, ynames, time_var, family,
                             age_var, age_center, control,
                             id_var = "id", rowid_var = "rowid",
                             maxit = 20L, tol = 1e-8, verbose = FALSE) {
  Q  <- length(ynames)
  Fa <- make_age_F(d[[age_var]] - age_center, d[[time_var]])
  ymat <- as.matrix(d[, ynames, drop = FALSE])

  gamma <- matrix(0, ncol(Fa), Q, dimnames = list(colnames(Fa), ynames))
  node  <- rep(NA_integer_, nrow(d))
  fml <- as.formula(paste0("cbind(", paste(ynames, collapse = ","), ") ~ ",
                           rowid_var, " + ", time_var, " | ",
                           paste(zvars, collapse = " + ")))

  tr <- NULL; conv <- FALSE
  for (it in seq_len(maxit)) {
    off  <- Fa %*% gamma
    fitf <- if (diagnostic == "slope")
      era_gee_mob_fit_slopeadj_off(d, ynames, time_var, family, off,
                                   rowid_var = rowid_var, id_var = id_var)
    else
      era_gee_mob_fit_off(d, ynames, time_var, family, off, rowid_var = rowid_var)

    tr <- mob(fml, data = d, fit = fitf, cluster = d[[id_var]], control = control)

    # membership mapped back through rowid; mob drops rows with NA on any Z
    mb <- rep(NA_integer_, nrow(d))
    mb[tr$data[[rowid_var]]] <- as.integer(predict(tr, type = "node"))
    ok <- !is.na(mb)

    g_new <- update_gamma(ymat[ok, , drop = FALSE], d[[time_var]][ok],
                          Fa[ok, , drop = FALSE], mb[ok])
    dg    <- max(abs(g_new - gamma))
    same  <- identical(mb, node)
    if (verbose)
      cat(sprintf("    iter %2d: %d terminal nodes, max|d gamma| = %.3g%s\n",
                  it, width(tr), dg, if (same) "  (membership stable)" else ""))
    gamma <- g_new; node <- mb
    if (same && dg < tol) { conv <- TRUE; break }
  }

  ok <- !is.na(node)
  list(tree = tr, gamma = gamma, node = node, iter = it, converged = conv,
       age_center = age_center,
       pernode_gamma = pernode_gamma(ymat[ok, , drop = FALSE],
                                     d[[time_var]][ok],
                                     Fa[ok, , drop = FALSE], node[ok]),
       offset = Fa %*% gamma)
}
