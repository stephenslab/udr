# Perform an M-step update for the mixture weights in the mixture
# prior.
update_mixture_weights_em <- function (P)
  colSums(P)/nrow(P)

# Perform an M-step update for the residual covariance matrix.
update_resid_covariance <- function (X, U, V, P) {
  n <- nrow(X)
  m <- ncol(X)
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
  Unew <- U
  if (version == "R") {
    k <- ncol(P)
    for (i in 1:k)
      Unew[,,i] <- update_prior_covariance_ed(X,U[,,i],V,P[,i])
  } else if (version == "Rcpp")
    Unew <- update_prior_covariances_ed_rcpp(X,U,V,P)
  return(Unew)
}

# Perform an M-step update for one of the prior covariance matrices
# using the update formula derived in Bovy et al (2011). Here, p is a
# vector, with one entry per row of X, giving the posterior assignment
# probabilities for the mixture component being updated.
update_prior_covariance_ed <- function (X, U, V, p) {
  p <- safenormalize(p)
  T <- U + V
  B <- solve(T,U)
  X1 <- crossprod((sqrt(p) * X) %*% B)
  return(U - U %*% B + X1)
}

# Perform an M-step update for the prior covariance matrices using the
# eigenvalue-truncation technique described in Won et al (2013). The
# calculations are implemented in both R (version = "R") and C++
# (version = "Rcpp").
update_prior_covariances_teem <- function (X, V, P, minval,
                                           version = c("Rcpp","R"), covtype) {
  version <- match.arg(version)
  if (version == "R") {
    m <- ncol(X)
    k <- ncol(P)
    U <- array(0,dim = c(m,m,k))
    for (i in 1:k)
      U[,,i] <- update_prior_covariance_teem(X,V,P[,i],minval, covtype)
  } else if (version == "Rcpp")
    U <- update_prior_covariances_teem_rcpp(X,V,P,minval)
  return(U)
}

# Perform an M-step update for one of the prior covariance matrices
# using the eigenvalue-truncation technique described in Won et al
# (2013). Input p is a vector, with one entry per row of X, giving
# the posterior assignment probabilities for the mixture components
# being updated.
update_prior_covariance_teem <- function (X, V, p, minval, covtype = c('unconstrained', 'rank1')) {

  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  p <- safenormalize(p)
  R <- chol(V)
  T <- crossprod((sqrt(p)*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite. This update for U is based on the fact
  # that the covariance matrix that minimizes the likelihood subject
  # to the constraint that U is positive definite is obtained by
  # truncating the eigenvalues of T = U + I less than 1 to be 1; see
  # Won et al (2013), p. 434, the sentence just after equation (16).
  if (covtype == 'unconstrained'){
        U <- shrink_cov(T, minval)
  }
  if (covtype == 'rank1'){
      U <- shrink_cov_deficient(T, r = 1, minval)
  }
  
  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}

# "Shrink" matrix T = U + I; that is, find the "best" matrix T
# satisfying the constraint that U is positive definite (if minval >
# 0) or positive semi-definite (if minval <= 0). This is achieved by
# setting any eigenvalues of T less than 1 to 1 + minval, or,
# equivalently, setting any eigenvalues of U less than 0 to be
# minval. The output is a positive definite matrix, U, or a positive
# semi-definite matrix if minval <= 0. 
shrink_cov <- function (T, minval = 0) {
  minval <- max(0,minval)
  out <- eigen(T)
  d <- out$values
  d <- pmax(d - 1,minval)
  return(tcrossprod(out$vectors %*% diag(sqrt(d))))
}


# r is the rank of U
shrink_cov_deficient <- function (T, r = 1, minval = 0) {
    
  minval <- max(0,minval)
  out <- eigen(T)
  d <- out$values
  d[(r+1):length(d)] <- 1   # keep only first r eigenvalues
  d <- pmax(d - 1,minval)
  return(tcrossprod(out$vectors %*% diag(sqrt(d))))
}


# function for 1-d search of s value based on eq.(20) in the write-up
# return a function of a scaler s_k.
# p: weights for component k
# s: scaling factor for k to search for
# Y: transformed data
# lambdas: eigenvalues of U_k
optimize_a_scaler = function(s, p, Y, lambdas){
  unweighted_sum = apply(Y, 2, function(y) sum(lambdas*y^2/((s*lambdas+1)^2))-sum(lambdas/(s*lambdas+1)))
  val = sum(p*unweighted_sum)
  return(val)
}
