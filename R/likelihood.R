# Compute the log-likelihood for the Ultimate Deconvolution model.
# Input argument U should be an m x m x k array, where m is the
# dimension of the data points, and k is the number of mixture
# components in the mixture-of-multivariate-normals prior. Input
# argument should either be an m x m matrix, or an m x m x n array,
# where n is the number of data points.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud <- function (X, w, U, V, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (is.list(U))
    U <- simplify2array(U)
  if (is.matrix(V)) {

    # Perform the computations for the case when the residual variance
    # is the same for all samples.
    if (version == "Rcpp")
      y <- loglik_ud_iid_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_iid_helper(X,w,U,V)
  } else {

    # Perform the computations for the case when the residual variance
    # is *not* the same for all samples.
    if (is.list(V))
      V <- simplify2array(V)
    if (version == "Rcpp")
      y <- loglik_ud_general_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_general_helper(X,w,U,V)
  }
  return(y)
}

# This is the R implementation of the likelihood computation when the
# residual covariance V is the same for all samples.
loglik_ud_iid_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (j in 1:k)
    y <- y + w[j] * dmvnorm(X,sigma = V + U[,,j])
  return(sum(log(y)))
}

# This is the R implementation of the likelihood computation when the
# residual covariance V is *not* the same for all samples.
loglik_ud_general_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (i in 1:n)
    for (j in 1:k)
      y[i] <- y[i] + w[j] * dmvnorm(X[i,],sigma = V[,,i] + U[,,j])
  return(sum(log(y)))
}
