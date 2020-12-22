# Perform an M-step update for the mixture weights in the
# mixture-of-multivariate-normals prior.
update_mixture_weights_em <- function (P)
  colSums(P)/nrow(P)

# Perform an M-step update for the residual covariance matrix.
update_resid_covariance <- function (X, U, V, P) {
  n    <- nrow(X)
  m    <- ncol(X)
  Vnew <- matrix(0,m,m)
  for (i in 1:n) {
    out  <- compute_posterior_mvtnorm_mix(X[i,],P[i,],U,V)
    Vnew <- Vnew + out$S1 + tcrossprod(X[i,] - out$mu1)
  }
  return(Vnew/n)
}

# Perform an M-step update for the prior covariance matrices using the
# update forumla derived in Bovy et al (2011). The calculations are
# implemented in both R (version = "R") and C++ (version = "Rcpp").
update_prior_covariances_ed <- function (X, U, V, P, version = c("Rcpp","R")) {
  version <- match.arg(version)
  U0 <- U
  if (version == "R") {
    k <- ncol(P)
    for (i in 1:k)
      U[,,i] <- update_prior_covariance_ed(X,U[,,i],V,P[,i])
  } else if (version == "Rcpp")
    U <- update_prior_covariances_ed_rcpp(X,U,V,P)
  return(U)
}

# Perform an M-step update for one of the prior covariance matrices
# using the update formula derived in Bovy et al (2011). Here, p is a
# vector, with one entry per row of X, giving the posterior assignment
# probabilities for the mixture component being updated.
update_prior_covariance_ed <- function (X, U, V, p) {
  T  <- U + V
  B  <- solve(T,U)
  X1 <- crossprod((sqrt(p/sum(p)) * X) %*% B)
  return(U - U %*% B + X1)
}

# Perform an M-step update for the prior covariance matrices using the
# eigenvalue-truncation technique described in Won et al (2013). The
# calculations are implemented in both R (version = "R") and C++
# (version = "Rcpp").
update_prior_covariances_teem <- function (X, V, P, minval,
                                           version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R") {
    m <- ncol(X)
    k <- ncol(P)
    U <- array(0,dim = c(m,m,k))
    for (i in 1:k)
      U[,,i] <- update_prior_covariance_teem(X,V,P[,i],minval)
  } else if (version == "Rcpp")
    U <- update_prior_covariances_teem_rcpp(X,V,P,minval)
  return(U)
}

# Perform an M-step update for one of the prior covariance matrices
# using the eigenvalue-truncation technique described in Won et al
# (2013). Input p is a vector, with one entry per row of X, giving
# the posterior assignment probabilities for the mixture components
# being updated.
update_prior_covariance_teem <- function (X, V, p, minval) {

  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  R <- chol(V)
  T <- crossprod((sqrt(p/sum(p))*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite. This update for U is based on the fact
  # that the covariance matrix that minimizes the likelihood subject
  # to the constraint that U is positive definite is obtained by
  # truncating the eigenvalues of T = U + I less than 1 to be 1; see
  # Won et al (2013), p. 434, the sentence just after equation (16).
  U <- shrink_cov(T,minval)

  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}
