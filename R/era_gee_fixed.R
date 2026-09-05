# ---------------------------------------------------------------------------
# era_gee_fixed.R
#
# ERA-GEE with FIXED components: the base model used throughout this project.
#
#     g(mu_itq) = beta_1q + beta_2q * t          q = 1, ..., Q
#
#     f_const = 1        component weight fixed at 1
#     f_time  = t        component weight fixed at 1, RAW wave index 0, 1, 2
#     B       = 2 x Q    row 1 = outcome-specific intercepts
#                        row 2 = outcome-specific time slopes
#
# Rows are person-waves. The GEE working covariance V_i is Q x Q: it models
# association among the Q outcomes WITHIN a person-wave. Wave-to-wave
# correlation is not parameterised; the fact that repeated observations belong
# to the same person is reflected instead in the cluster-level estimating
# functions and in robust instability inference, via mob(..., cluster = id).
#
# ---------------------------------------------------------------------------
# WHY THE COMPONENTS ARE FIXED
#
# Under this specification F = [1  t] is a GIVEN design matrix, not something
# to be estimated. Two things therefore have to be switched off relative to a
# general ERA-GEE implementation:
#
#   * the weight update. w_const = w_time = 1 are CONSTRAINTS. "Converges to
#     approximately 1" and "fixed at 1" are different models.
#   * the normalisation. Standardising X would make a constant column
#     identically zero, destroying the intercept component; standardising F
#     would turn t into (t - tbar)/s_t, after which beta_2q is no longer
#     "change per wave".
#
# What remains is a multivariate GEE with a fixed design. `ERA_GEE_fixed()`
# asserts the invariants max|f_const - 1| = 0 and max|f_time - t| = 0 on every
# call.
#
# PARAMETER ORDER. vec(t(B)) = (beta_11 ... beta_1Q, beta_21 ... beta_2Q):
# all intercepts, then all slopes. So estfun columns 1:Q are the intercept
# scores and (Q+1):(2Q) the slope scores.
#
# OUTCOME NORMALISATION. `y_scale` is applied once and NEVER per node:
#   "raw"    outcomes on their original scale (used throughout this project;
#            intercepts and slopes are directly interpretable)
#   "global" centre/scale once on the full analytic sample, holding those
#            constants fixed at every node (pass y_center and y_sd)
# Node-specific re-standardisation is not offered: it is what breaks the
# comparability of beta across terminal nodes.
#
# ---------------------------------------------------------------------------
# SCOPE OF THIS IMPLEMENTATION
#
# Gaussian outcomes with the identity link, and `corstr = "independence"`.
# That is the specification that was validated (see docs/VALIDATION.md); other
# families and working correlation structures are rejected with an error
# rather than silently approximated.
#
# Dependencies: MASS (ginv, mvrnorm), matrixcalc (vec). No other code needed.
# ---------------------------------------------------------------------------


# --- family step ----------------------------------------------------------
# Adjusted dependent variable, mean and variance function.
# For the Gaussian identity link: mu = eta, mvar = 1, and the adjusted
# dependent variable z = eta + (y - mu)/(dmu/deta) = y.
gee_family_step <- function(y, eta, family = "gaussian") {
  if (!identical(family, "gaussian"))
    stop("era_gee_fixed.R implements the Gaussian identity link only; got '",
         family, "'. See docs/VALIDATION.md for the validated scope.")
  list(z = as.numeric(y), mu = as.numeric(eta), mvar = rep(1, length(eta)))
}

# --- working precision ----------------------------------------------------
# Inverse working covariance V_i^{-1} for one person-wave cluster of Q
# outcomes. With `corstr = "independence"` the working correlation is I_Q, and
# for the Gaussian identity link the variance function is 1, so
# V_i = A^{1/2} R A^{1/2} = I_Q and V_i^{-1} = I_Q for every i.
#
# The dispersion phi is returned for reporting but does not enter V_i^{-1}:
# it cancels out of both the estimating equation and the score, so including
# it would rescale the objective without changing any estimate or test.
gee_working_precision <- function(resid, mvar, Q, corstr = "independence") {
  if (!identical(corstr, "independence"))
    stop("era_gee_fixed.R implements corstr = 'independence' only; got '",
         corstr, "'. See docs/VALIDATION.md for the validated scope.")
  phi <- sum(resid^2) / length(resid)
  list(Vinv = diag(Q), R = diag(Q), phi = phi)
}


# --- main fitter ----------------------------------------------------------
ERA_GEE_fixed <- function(y, time, family, corstr = "independence",
                          y_scale = c("raw", "global"),
                          y_center = NULL, y_sd = NULL,
                          maxit = 100, ceps = 1e-5, verbose = FALSE) {

  y_scale <- match.arg(y_scale)
  y <- as.matrix(y)
  ynames <- colnames(y)
  if (is.null(ynames)) ynames <- paste0("Y", seq_len(ncol(y)))

  nobs <- nrow(y)      # person-wave rows
  Q    <- ncol(y)      # outcomes -> GEE cluster size
  nlv  <- 2L           # intercept component + time component

  if (length(family) != Q) stop("length(family) must equal ncol(y)")
  if (length(time)   != nobs) stop("length(time) must equal nrow(y)")

  if (y_scale == "global") {
    if (is.null(y_center) || is.null(y_sd))
      stop('y_scale = "global" requires y_center and y_sd from the full sample')
    y <- sweep(sweep(y, 2, y_center, "-"), 2, y_sd, "/")
  }

  # FIXED components. Never rescaled, never re-estimated.
  F <- cbind(f_const = rep(1, nobs), f_time = as.numeric(time))
  stopifnot(max(abs(F[, "f_const"] - 1))    == 0,
            max(abs(F[, "f_time"]  - time)) == 0)

  B <- MASS::ginv(t(F) %*% F) %*% t(F) %*% y            # OLS start
  dimnames(B) <- list(c("intercept", "slope"), ynames)

  Ivec <- diag(Q)
  conv <- FALSE

  for (iter in seq_len(maxit)) {
    lp <- F %*% B
    Z <- mu <- mvar <- matrix(NA_real_, nobs, Q)
    for (r in seq_len(Q)) {
      st <- gee_family_step(y[, r], lp[, r], family[r])
      Z[, r] <- st$z; mu[, r] <- st$mu; mvar[, r] <- st$mvar
    }
    wp <- gee_working_precision(as.numeric(y - mu), mvar, Q, corstr)

    # B update: weighted least squares with the FIXED design
    T1 <- matrix(0, nlv * Q, nlv * Q); T2 <- matrix(0, nlv * Q, 1)
    for (i in seq_len(nobs)) {
      Mi  <- kronecker(F[i, , drop = FALSE], Ivec)      # Q x 2Q
      tMV <- t(Mi) %*% wp$Vinv
      T1  <- T1 + tMV %*% Mi
      T2  <- T2 + tMV %*% matrix(Z[i, ], ncol = 1)
    }
    B_new <- t(matrix(MASS::ginv(T1) %*% T2, nrow = Q, ncol = nlv))
    dimnames(B_new) <- list(c("intercept", "slope"), ynames)

    delta <- sum(abs(B_new - B)); B <- B_new
    if (delta < ceps) { conv <- TRUE; break }
  }
  if (verbose)
    cat(if (conv) sprintf("converged in %d iterations\n", iter)
        else sprintf("NOT converged in %d iterations\n", maxit))

  # final pass: objective function and score contributions
  lp <- F %*% B
  Z <- mu <- mvar <- matrix(NA_real_, nobs, Q)
  for (r in seq_len(Q)) {
    st <- gee_family_step(y[, r], lp[, r], family[r])
    Z[, r] <- st$z; mu[, r] <- st$mu; mvar[, r] <- st$mvar
  }
  wp <- gee_working_precision(as.numeric(y - mu), mvar, Q, corstr)

  vecB   <- as.numeric(matrixcalc::vec(t(B)))           # intercepts then slopes
  objfun <- 0
  estfun <- matrix(NA_real_, nobs, nlv * Q)
  for (i in seq_len(nobs)) {
    fi <- F[i, , drop = FALSE]
    ri <- matrix(Z[i, ] - as.numeric(fi %*% B), ncol = 1)          # Q x 1
    objfun      <- objfun + as.numeric(t(ri) %*% wp$Vinv %*% ri)
    estfun[i, ] <- -2 * t(kronecker(fi, Ivec)) %*% wp$Vinv %*% ri
  }
  colnames(estfun) <- c(paste0("b1_", ynames), paste0("b2_", ynames))

  list(B = B, coef = vecB, objfun = objfun, estfun = estfun,
       F = F, phi = wp$phi, R = wp$R, corstr = corstr,
       converged = conv, iter = iter, y_scale = y_scale,
       ynames = ynames, Q = Q)
}


# ---------------------------------------------------------------------------
# mob() fit-function factory -- METHOD A (the one used for all results).
#
# Node membership is read from x[, "rowid"], NOT from `weights`. partykit
# passes weights already subset to the node (all positive) at every node below
# the root, so `which(weights > 0)` returns 1:n_node and would silently select
# the first n_node rows of the whole dataset.
#
# Requires `rowid` on the regressor side of the mob formula:
#     cbind(y1,...,yQ) ~ rowid + X_time | z1 + z2 + ...
# and mob(..., cluster = id) for person-level instability inference.
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
    if (is.null(x)) stop("needs `x`; put `rowid` on the regressor side")
    xm <- as.matrix(x)
    if (!rowid_var %in% colnames(xm))
      stop("`", rowid_var, "` not in the model matrix; add it to the formula")

    d <- dat.full[as.integer(round(xm[, rowid_var])), , drop = FALSE]
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
# mob() fit-function factory -- METHOD B(a): subject-level interface.
#
# The theoretical reference implementation. In GEE the independent unit is the
# subject and the natural estimating-function contribution is
# U_i = sum_t u_it. Because the working covariance is block-diagonal across
# waves, that decomposition is exact and unique. Here mob() partitions a
# PERSON-level frame, so estfun is N x 2Q and the instability tests run over
# genuinely independent clusters. No `cluster` argument is needed.
#
# Used in run/02_cluster_equivalence.R to check that method A reproduces it.
#
# `pdat` needs `prow` (1..N) on the regressor side and a column `id_var`
# matching `dat.full`. Its Q outcome columns are placeholders for mob()'s
# model frame and are never read.
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
    node_ids <- ids_by_prow[as.integer(round(xm[, prow_var]))]

    d <- dat.full[unlist(split_rows[as.character(node_ids)], use.names = FALSE), ,
                  drop = FALSE]
    m <- ERA_GEE_fixed(y = d[, ynames, drop = FALSE], time = d[[time_var]],
                       family = family, corstr = corstr,
                       y_scale = y_scale, y_center = y_center, y_sd = y_sd)

    rval <- list(coefficients = setNames(m$coef, colnames(m$estfun)),
                 objfun       = m$objfun)
    if (estfun) {
      ef <- rowsum(m$estfun, group = d[[id_var]], reorder = FALSE)   # U_i = sum_t u_it
      ef <- ef[match(as.character(node_ids), rownames(ef)), , drop = FALSE]
      rownames(ef) <- NULL
      rval$estfun <- ef
    }
    if (object) rval$object <- m
    rval
  }
}
