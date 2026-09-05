# ---------------------------------------------------------------------------
# nuisance_adjusted_slope.R
#
# PROTOTYPE: nuisance-orthogonalized slope score for the slope-focused
# diagnostic. Existing factories in ERA_GEE_fixed.R are left untouched.
#
# CONCEPTUAL SPECIFICATION (this much is fixed):
#   Fit the FULL 2Q model. Partition the estimating-function contribution as
#   U_i = (U_{I,i}, U_{S,i}), intercept block and slope block. Perform the
#   instability test on the Q nuisance-adjusted slope scores
#
#       U_{S.I} = U_S - I_{SI} I_{II}^{-1} U_I
#
#   i.e. project the nuisance intercept direction out of the slope score.
#
# Everything else -- what `coefficients` and `objfun` should carry, what the
# terminal-node return structure should look like -- is an interface question
# for this prototype to answer, NOT a theoretical design decision. The values
# reported below are placeholders chosen to satisfy mob()'s interface.
#
# NAMING. "nuisance-adjusted" / "nuisance-orthogonalized", deliberately not
# "efficient score": GEE is built on estimating equations rather than a full
# likelihood, so calling this efficient would need extra justification about
# which empirical covariance / sensitivity matrix is being used. The
# time-origin invariance below is a purely algebraic property and can be
# tested as it stands.
#
# WHY THIS SHOULD BE INVARIANT TO THE TIME ORIGIN.
#   Under t -> t - c the design transforms as F_new = F A with
#   A = [[1, -c], [0, 1]], so the scores transform as (A' kron I_Q):
#
#       U_{I,new} = U_I ,      U_{S,new} = U_S - c U_I .
#
#   The slope score changes only by a multiple of the intercept score, which is
#   exactly the direction the projection removes:
#
#       resid(U_S - c U_I | U_I) = resid(U_S | U_I) .
#
#   So U_{S.I} is invariant elementwise, and the resulting instability
#   statistic must agree under raw and centred time to numerical tolerance.
#   PASS CRITERION for this prototype. A failure indicts the implementation
#   first, not the method.
#
# WHERE THE PROJECTION IS ESTIMATED.
#   I_{SI} I_{II}^{-1} is estimated from the SUBJECT-level scores U_i within
#   the current node (the clustered test's unit), then applied row-wise to the
#   wave-level scores. Because the projection is linear, summing the adjusted
#   wave-level scores over a person recovers the adjusted subject-level score,
#   so method A + cluster = id stays coherent with B(a).
#
# NOTE: `parm` is NOT used. partykit whitens the process and then subsets its
# columns, which would re-mix the blocks; here the Q-dimensional adjusted
# process is handed over directly and whitened on its own covariance.
# ---------------------------------------------------------------------------

# Returns the wave-level nuisance-adjusted slope scores (n_row x Q) for a
# fitted ERA_GEE_fixed object, given the person id of each row.
nuisance_adjust_slope <- function(estfun, id, Q) {
  UI <- estfun[, 1:Q, drop = FALSE]
  US <- estfun[, (Q + 1):(2 * Q), drop = FALSE]

  # projection estimated on the SUBJECT-level scores
  A  <- rowsum(estfun, group = id, reorder = FALSE)
  AI <- A[, 1:Q, drop = FALSE]
  AS <- A[, (Q + 1):(2 * Q), drop = FALSE]
  Bmat <- MASS::ginv(crossprod(AI)) %*% crossprod(AI, AS)      # Q x Q

  adj <- US - UI %*% Bmat
  colnames(adj) <- paste0("b2adj_", sub("^b1_", "", colnames(UI)))
  adj
}

era_gee_mob_fit_slopeadj <- function(dat.full, ynames, time_var, family,
                                     corstr = "independence",
                                     y_scale = c("raw", "global"),
                                     y_center = NULL, y_sd = NULL,
                                     rowid_var = "rowid", id_var = "id") {
  y_scale <- match.arg(y_scale)
  Q <- length(ynames)

  function(y, x = NULL, start = NULL, weights = NULL, offset = NULL,
           estfun = FALSE, object = FALSE, ...) {
    xm <- as.matrix(x)
    node_rows <- as.integer(round(xm[, rowid_var]))
    d <- dat.full[node_rows, , drop = FALSE]

    # the FULL 2Q model is what is actually fitted
    m <- ERA_GEE_fixed(y = d[, ynames, drop = FALSE], time = d[[time_var]],
                       family = family, corstr = corstr,
                       y_scale = y_scale, y_center = y_center, y_sd = y_sd)

    adj <- nuisance_adjust_slope(m$estfun, d[[id_var]], Q)

    # INTERFACE EXPEDIENT, not a design decision: mob() wants `coefficients`
    # to line up with the columns of `estfun`, so the slope half is reported
    # here. The full 2 x Q model is still what was fitted, and is available
    # via `object`.
    rval <- list(coefficients = setNames(m$B["slope", ], colnames(adj)),
                 objfun       = m$objfun)
    if (estfun) rval$estfun <- adj
    if (object) rval$object <- m
    rval
  }
}
