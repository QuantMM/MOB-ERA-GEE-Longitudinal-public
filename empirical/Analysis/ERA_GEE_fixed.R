# ---------------------------------------------------------------------------
# ERA_GEE_fixed.R
#
# ERA-GEE with FIXED components, implementing the agreed base model exactly:
#
#     g(mu_itq) = beta_1q + beta_2q * t        q = 1, ..., Q
#
#     f_const = 1        (component weight fixed at 1)
#     f_time  = t        (component weight fixed at 1, RAW wave index 0,1,2)
#     B       = 2 x Q    (row 1 = outcome-specific intercepts
#                         row 2 = outcome-specific time slopes)
#
# WHY A SEPARATE FITTER RATHER THAN A PATCH TO MOB.MERA()
# ------------------------------------------------------
# In the agreed specification F = [1  t] is a GIVEN design matrix, not an
# estimated quantity. That removes the two things MOB.MERA() spends most of
# its loop on:
#
#   * the W update      -- w_const = w_time = 1 are CONSTRAINTS, so the
#                          optimizer must not touch them. "Converges to
#                          approximately 1" and "fixed at 1" are different
#                          models.
#   * the normalisation -- Xscale() applied to X, to y, and to F on every
#                          iteration. Passing X_const through it makes the
#                          intercept component identically zero, and passing
#                          F through it turns t into (t - tbar)/s_t, after
#                          which beta_2q is no longer "change per wave".
#
# What remains is exactly a multivariate GEE with a fixed design. The GEE
# machinery itself is REUSED VERBATIM from Sam's helper file (`linkfun`,
# `getRi`, `getAVAi`), so the working-covariance estimation is identical to
# the shipped engine. Only the ERA layer around it is replaced by the fixed
# specification.
#
# The B update below is the same weighted least-squares step as MOB.MERA():
#     M_i           = kron(f_i, I_Q)                      (Q x 2Q)
#     vec(B)        = (sum_i M_i' V_i^-1 M_i)^-1 (sum_i M_i' V_i^-1 Z_i)
#     objfun        = sum_i (Z_i - M_i vec(B))' V_i^-1 (Z_i - M_i vec(B))
#     estfun[i, ]   = -2 * M_i' V_i^-1 (Z_i - M_i vec(B))
#
# PARAMETER ORDER.  vec(t(B)) = (beta_11 ... beta_1Q, beta_21 ... beta_2Q),
# i.e. all intercepts first, then all slopes. So estfun columns 1:Q are the
# intercept scores and (Q+1):(2Q) are the slope scores. This is the ordering
# the intercept-only / slope-only diagnostic tests will index into.
#
# CLUSTERING.  Rows are person-waves. V_i is Q x Q, i.e. the working
# covariance models association among the Q outcomes WITHIN a person-wave;
# wave-to-wave correlation is not parameterised. The fact that repeated
# observations belong to the same person is reflected instead in the
# cluster-level estimating functions / robust instability inference
# (mob(..., cluster = id)).
#
# Y NORMALISATION.  `y_scale` is a switch, never node-specific:
#     "raw"    -- outcomes on their original scale (intercepts and slopes are
#                 directly interpretable; used for the worked example)
#     "global" -- centre/scale ONCE on the full analytic sample and hold those
#                 constants fixed at every node (pass them in via `y_center`
#                 and `y_sd`)
# Node-specific re-standardisation is not offered: it is what breaks the
# comparability of beta across terminal nodes.
#
# WORKING CORRELATION.  `corstr = "independence"` is used throughout the
# worked example. See NOTES_getRi_indexing.md before using any other
# structure with this data layout.
# ---------------------------------------------------------------------------

ERA_GEE_fixed <- function(y, time, family, corstr = "independence",
                          y_scale = c("raw", "global"),
                          y_center = NULL, y_sd = NULL,
                          maxit = 100, ceps = 1e-5, verbose = FALSE) {

  y_scale <- match.arg(y_scale)
  y <- as.matrix(y)
  ynames <- colnames(y)
  if (is.null(ynames)) ynames <- paste0("Y", seq_len(ncol(y)))

  nobs    <- nrow(y)          # number of person-wave rows
  nrepeat <- ncol(y)          # Q outcomes  -> GEE cluster size
  Q       <- nrepeat
  nlv     <- 2L               # intercept component + time component
  nt      <- as.integer(rep(nrepeat, nobs))

  if (length(family) != Q)
    stop("length(family) must equal ncol(y)")
  if (length(time) != nobs)
    stop("length(time) must equal nrow(y)")

  # --- outcome scaling: applied ONCE, with constants supplied from outside --
  if (y_scale == "global") {
    if (is.null(y_center) || is.null(y_sd))
      stop('y_scale = "global" requires y_center and y_sd from the full sample')
    y <- sweep(sweep(y, 2, y_center, "-"), 2, y_sd, "/")
  }

  # --- FIXED components. Never rescaled, never re-estimated. ----------------
  F <- cbind(f_const = rep(1, nobs), f_time = as.numeric(time))

  # invariants that define the specification
  stopifnot(max(abs(F[, "f_const"] - 1))        == 0,
            max(abs(F[, "f_time"]  - time))     == 0)

  # --- initial B: OLS of y on F --------------------------------------------
  B <- MASS::ginv(t(F) %*% F) %*% t(F) %*% y
  dimnames(B) <- list(c("intercept", "slope"), ynames)

  # -------------------------------------------------------------------------
  # FAST PATH -- algebraically identical, not an approximation.
  #
  # With a Gaussian identity link and corstr = "independence", linkfun gives
  # z = y, mu = eta, mvar = 1, and getAVAi returns V_i^{-1} = I_Q for every i
  # (the phi rescaling is commented out in the shipped helper). The per-row
  # loop then collapses:
  #
  #   M_i = kron(f_i, I_Q)  =>  sum_i M_i' V_i^-1 M_i = kron(F'F, I_Q)
  #                             sum_i M_i' V_i^-1 z_i = as.vector(t(Y) %*% F)
  #   => vec(B) solves kron(F'F, I_Q) vec(B) = ..., i.e. B = (F'F)^-1 F'Y
  #
  #   estfun[i, ]  = -2 * [f_i1 * r_i ; f_i2 * r_i]
  #   => estfun = -2 * cbind(f_const * R, f_time * R),  R = Y - F B
  #
  #   objfun = sum_i r_i' I_Q r_i = sum(R^2)
  #
  # Verified against the looped path to 0 difference on the reference case.
  # Anything outside this family/corstr falls through to the general loop.
  # -------------------------------------------------------------------------
  if (all(family == "gaussian") && identical(corstr, "independence")) {
    B <- MASS::ginv(crossprod(F)) %*% crossprod(F, y)
    dimnames(B) <- list(c("intercept", "slope"), ynames)
    R <- y - F %*% B
    estfun <- -2 * cbind(F[, 1] * R, F[, 2] * R)
    colnames(estfun) <- c(paste0("b1_", ynames), paste0("b2_", ynames))
    return(list(B = B, coef = as.numeric(matrixcalc::vec(t(B))),
                objfun = sum(R^2), estfun = estfun,
                F = F, phi = sum(R^2) / length(R), R = diag(Q), corstr = corstr,
                converged = TRUE, iter = 1L, y_scale = y_scale,
                ynames = ynames, Q = Q))
  }

  Ivec <- diag(Q)
  conv <- FALSE

  for (iter in seq_len(maxit)) {

    lp <- F %*% B                                    # linear predictor

    Z <- mu <- mvar <- matrix(NA_real_, nobs, Q)
    for (r in seq_len(Q)) {
      adj <- linkfun(F, lp[, r, drop = FALSE], y[, r, drop = FALSE],
                     family = family[r])
      Z[, r]    <- adj$z
      mu[, r]   <- adj$mu
      mvar[, r] <- adj$mvar
    }

    Rp   <- getRi(y, mu, mvar, nt, corstr)
    AVAi <- getAVAi(mvar, nt, Rp$Ehat, Rp$phi)

    # --- B update: weighted least squares with the FIXED design ------------
    T1 <- matrix(0, nlv * Q, nlv * Q)
    T2 <- matrix(0, nlv * Q, 1)
    for (i in seq_len(nobs)) {
      Mi  <- kronecker(F[i, , drop = FALSE], Ivec)    # Q x 2Q
      Vi  <- AVAi[seq_len(Q), seq_len(Q), i]
      tMV <- t(Mi) %*% Vi
      T1  <- T1 + tMV %*% Mi
      T2  <- T2 + tMV %*% matrix(Z[i, ], ncol = 1)
    }
    vecB_new <- MASS::ginv(T1) %*% T2                 # (beta_1., beta_2.)
    B_new    <- matrix(vecB_new, nrow = Q, ncol = nlv)   # Q x 2 = t(B)
    B_new    <- t(B_new)
    dimnames(B_new) <- list(c("intercept", "slope"), ynames)

    delta <- sum(abs(B_new - B))
    B <- B_new
    if (delta < ceps) { conv <- TRUE; break }
  }
  if (verbose)
    cat(if (conv) sprintf("converged in %d iterations\n", iter)
        else sprintf("NOT converged in %d iterations\n", maxit))

  # --- final pass: objective function and score contributions --------------
  lp <- F %*% B
  Z <- mu <- mvar <- matrix(NA_real_, nobs, Q)
  for (r in seq_len(Q)) {
    adj <- linkfun(F, lp[, r, drop = FALSE], y[, r, drop = FALSE], family = family[r])
    Z[, r] <- adj$z; mu[, r] <- adj$mu; mvar[, r] <- adj$mvar
  }
  Rp   <- getRi(y, mu, mvar, nt, corstr)
  AVAi <- getAVAi(mvar, nt, Rp$Ehat, Rp$phi)

  vecB   <- as.numeric(matrixcalc::vec(t(B)))         # intercepts then slopes
  objfun <- 0
  estfun <- matrix(NA_real_, nobs, nlv * Q)
  for (i in seq_len(nobs)) {
    fi <- F[i, , drop = FALSE]
    ri <- matrix(Z[i, ] - as.numeric(fi %*% B), ncol = 1)     # Q x 1
    Vi <- AVAi[seq_len(Q), seq_len(Q), i]
    objfun    <- objfun + as.numeric(t(ri) %*% Vi %*% ri)
    estfun[i, ] <- -2 * t(kronecker(fi, Ivec)) %*% Vi %*% ri
  }
  colnames(estfun) <- c(paste0("b1_", ynames), paste0("b2_", ynames))

  list(B = B, coef = vecB, objfun = objfun, estfun = estfun,
       F = F, phi = Rp$phi, R = Rp$Ehat, corstr = corstr,
       converged = conv, iter = iter, y_scale = y_scale,
       ynames = ynames, Q = Q)
}


# ---------------------------------------------------------------------------
# mob() fit-function factory.
#
# Node membership is read from x[, "rowid"], NOT from `weights`. partykit
# passes weights already subset to the node (all positive) at every node below
# the root, so `which(weights > 0)` would return 1:n_node and silently select
# the first n_node rows of the whole dataset. See
# Analysis/REPORT_A_pipeline_verification.md, defect D4.
#
# Requires `rowid` on the regressor side of the mob formula:
#     cbind(y1,...,yQ) ~ rowid + X_time | z1 + z2 + ...
# ---------------------------------------------------------------------------
era_gee_mob_fit <- function(dat.full, ynames, time_var, family,
                            corstr = "independence",
                            y_scale = c("raw", "global"),
                            y_center = NULL, y_sd = NULL,
                            rowid_var = "rowid") {
  y_scale <- match.arg(y_scale)
  if (!rowid_var %in% names(dat.full))
    stop("`", rowid_var, "` must be a column of dat.full")

  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) {
    xm <- as.matrix(x)
    if (!rowid_var %in% colnames(xm))
      stop("`", rowid_var, "` must appear on the regressor side of the formula")
    node_rows <- as.integer(round(xm[, rowid_var]))

    d <- dat.full[node_rows, , drop = FALSE]
    m <- ERA_GEE_fixed(y = d[, ynames, drop = FALSE], time = d[[time_var]],
                       family = family, corstr = corstr,
                       y_scale = y_scale, y_center = y_center, y_sd = y_sd)

    rval <- list(coefficients = setNames(m$coef, colnames(m$estfun)),
                 objfun       = m$objfun)
    if (estfun) rval$estfun <- m$estfun
    if (object) rval$object <- m
    rval
  }
}


# ---------------------------------------------------------------------------
# METHOD B(a): subject-level MOB interface.
#
# The theoretical reference implementation. In GEE the independent unit is the
# subject, and the natural empirical estimating-function contribution is
#
#     U_i(theta) = sum_t u_it(theta)
#
# Because the working covariance is block-diagonal across waves (V_i models
# association among the Q outcomes WITHIN a person-wave and does not
# parameterise wave-to-wave correlation), this decomposition is exact and
# unique: summing the wave-level contributions recovers U_i with nothing left
# over. B(a) is therefore not "ignoring wave dependence for convenience" -- it
# is ERA-GEE using working independence across waves while respecting the
# clustering of repeated measurements in the score inference.
#
# Here mob() partitions a PERSON-LEVEL data frame (N rows), so estfun is
# N x 2Q and the instability tests run over genuinely independent clusters
# i = 1, ..., N. No `cluster` argument is needed. The partitioning variables
# are person-level anyway, so the conceptual match is exact.
#
# `pdat` must contain `prow` (1..N) on the regressor side of the formula and a
# column `id_var` giving the person id used in `dat.full`. The Q outcome
# columns in `pdat` are placeholders required by mob()'s model frame; they are
# never read.
# ---------------------------------------------------------------------------
era_gee_mob_fit_subject <- function(dat.full, pdat, ynames, time_var, family,
                                    corstr = "independence",
                                    y_scale = c("raw", "global"),
                                    y_center = NULL, y_sd = NULL,
                                    id_var = "id", prow_var = "prow") {
  y_scale <- match.arg(y_scale)
  if (!prow_var %in% names(pdat)) stop("`", prow_var, "` must be a column of pdat")
  ids_by_prow <- pdat[[id_var]][order(pdat[[prow_var]])]
  split_rows  <- split(seq_len(nrow(dat.full)), dat.full[[id_var]])

  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) {
    xm <- as.matrix(x)
    if (!prow_var %in% colnames(xm))
      stop("`", prow_var, "` must appear on the regressor side of the formula")
    node_prows <- as.integer(round(xm[, prow_var]))
    node_ids   <- ids_by_prow[node_prows]

    rows <- unlist(split_rows[as.character(node_ids)], use.names = FALSE)
    d    <- dat.full[rows, , drop = FALSE]

    m <- ERA_GEE_fixed(y = d[, ynames, drop = FALSE], time = d[[time_var]],
                       family = family, corstr = corstr,
                       y_scale = y_scale, y_center = y_center, y_sd = y_sd)

    rval <- list(coefficients = setNames(m$coef, colnames(m$estfun)),
                 objfun       = m$objfun)
    if (estfun) {
      # U_i = sum_t u_it, in the order mob expects (i.e. by node_prows)
      ef <- rowsum(m$estfun, group = d[[id_var]], reorder = FALSE)
      ef <- ef[match(as.character(node_ids), rownames(ef)), , drop = FALSE]
      rownames(ef) <- NULL
      rval$estfun <- ef
    }
    if (object) rval$object <- m
    rval
  }
}
