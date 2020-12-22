# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance V, and z is drawn from a mixture of multivariate
# normals, each with zero mean, covariance U[,,i] and weight w[i].
# Return the posterior mean (mu1) and covariance (S1) of z. Note that
# input w1 must be the vector of *posterior* mixture weights (see
# compute_posterior_probs).
compute_posterior_mvtnorm_mix <- function (x, w1, U, V) {
  m   <- length(x)
  k   <- length(w1)
  mu1 <- rep(0,m)
  S1  <- matrix(0,m,m)
  for (i in 1:k) {
    out <- compute_posterior_mvtnorm(x,U[,,i],V)
    mu1 <- mu1 + w1[i] * out$mu1
    S1  <- S1 + w1[i] * (out$S1 + tcrossprod(out$mu1))
  }
  S1 <- S1 - tcrossprod(mu1)
  return(list(mu1 = mu1,S1 = S1))
}

# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance V, and z is drawn from a multivariate normal
# distribution with mean zero and covariance U. Return the posterior
# mean (mu) and covariance (S1) of z. These calculations will only
# work if V is positive definite (invertible).
compute_posterior_mvtnorm <- function (x, U, V) {
  m   <- length(x)
  S1  <- solve(U %*% solve(V) + diag(m)) %*% U
  mu1 <- drop(S1 %*% solve(V,x))
  return(list(mu1 = mu1,S1 = S1))
}

# This implements the calculations for compute_posterior_probs for the
# special case when the same residual covariance matrix is used for
# all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_iid <- function (X, w, U, V) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  k <- length(w)

  # Compute the log-probabilities, stored in an n x k matrix.
  P <- matrix(0,n,k)
  for (j in 1:k)
    P[,j] = log(w[j]) + dmvnorm(X,sigma = V + U[,,j],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(softmax(P))
}

# This implements the calculations for compute_posterior_probs for the
# more general case when the residual covariances are not necessaily
# the same for all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_general <- function (X, w, U, V) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  k <- length(w)

  # Compute the log-probabilities, stored in an n x k matrix.
  P <- matrix(0,n,k)
  for (i in 1:n)
    for (j in 1:k)
      P[i,j] = log(w[j]) + dmvnorm(X[i,],sigma = V[,,i] + U[,,j],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(softmax(P))
}

# Compute the n x k matrix of posterior mixture assignment
# probabilities given current estimates of the model parameters. This
# implements the E step in the EM algorithm for fitting the Ultimate
# Deconvolution model.
compute_posterior_probs <- function (X, w, U, V, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (is.matrix(V)) {

    # Perform the computations for the special case when the same
    # residual variance is used for all samples.
    if (version == "R")
      P <- compute_posterior_probs_iid(X,w,U,V)
    else if (version == "Rcpp")
      P <- compute_posterior_probs_iid_rcpp(X,w,U,V)
  } else {

    # Perform the computations for the more general case when the
    # residual variance is not necessarily the same for all samples.
    if (version == "R")
      P <- compute_posterior_probs_general(X,w,U,V)
    else if (version == "Rcpp")
      P <- compute_posterior_probs_general_rcpp(X,w,U,V)
  }
  return(P)
}
