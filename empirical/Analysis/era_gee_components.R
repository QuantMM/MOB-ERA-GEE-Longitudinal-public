# ---------------------------------------------------------------------------
# era_gee_components.R
#
# ERA-GEE with an arbitrary set of FIXED components, for the Gaussian family
# with an identity link and working independence:
#
#     g(mu_itq) = sum_k  b_kq * f_it,k        k = 1..K,  q = 1..Q
#
# The design F = [f_1 ... f_K] is given (e.g. [1, t], [1, t, A_c],
# [1, t, A_c, t*A_c]); every component weight is fixed at 1 and nothing is
# rescaled. B is K x Q and theta = vec(t(B)) orders all coefficients of
# component 1 first, then component 2, and so on -- the same convention as
# ERA_GEE_fixed(), so with F = [1, t] the two functions agree exactly.
#
# ESTIMATION. With identity link and working independence the GEE estimating
# equations are the normal equations, so B = (F'F)^-1 F'Y. Point estimates are
# therefore those of Q outcome-wise least-squares fits; what the joint
# formulation adds is the covariance of theta ACROSS outcomes, which is what
# multi-outcome Wald tests need.
#
# INFERENCE. Person-clustered sandwich (the usual GEE robust variance, no
# small-sample correction):
#     bread = (F'F)^-1 kron I_Q
#     meat  = sum_i U_i U_i',   U_i = sum_t  f_it kron r_it
#     V     = bread meat bread
# Repeated waves of the same person enter only through U_i, i.e. working
# independence across waves with person-level robust inference.
# ---------------------------------------------------------------------------

era_gee_components <- function(Y, F, id) {
  Y <- as.matrix(Y); F <- as.matrix(F)
  Q <- ncol(Y); K <- ncol(F)
  if (is.null(colnames(F))) colnames(F) <- paste0("f", seq_len(K))
  if (is.null(colnames(Y))) colnames(Y) <- paste0("Y", seq_len(Q))
  stopifnot(nrow(F) == nrow(Y), length(id) == nrow(Y), !anyNA(Y), !anyNA(F))

  FtF_inv <- solve(crossprod(F))
  B <- FtF_inv %*% crossprod(F, Y)
  dimnames(B) <- list(colnames(F), colnames(Y))
  R <- Y - F %*% B

  U  <- do.call(cbind, lapply(seq_len(K), function(k) F[, k] * R))   # n x KQ
  Ui <- rowsum(U, group = id, reorder = FALSE)                       # persons x KQ
  bread <- kronecker(FtF_inv, diag(Q))
  V <- bread %*% crossprod(Ui) %*% bread
  nm <- as.vector(t(outer(colnames(F), colnames(Y), paste, sep = ":")))
  dimnames(V) <- list(nm, nm)

  SE <- matrix(sqrt(diag(V)), K, Q, byrow = TRUE, dimnames = dimnames(B))
  list(B = B, V = V, SE = SE, theta = as.vector(t(B)), R = R,
       n_rows = nrow(Y), n_persons = nrow(Ui), K = K, Q = Q)
}

# Wald test that a set of elements of theta (named "component:outcome") is 0
era_wald <- function(fit, which) {
  b <- fit$theta[match(which, rownames(fit$V))]
  W <- as.numeric(t(b) %*% solve(fit$V[which, which, drop = FALSE]) %*% b)
  c(W = W, df = length(which), p = pchisq(W, length(which), lower.tail = FALSE))
}
