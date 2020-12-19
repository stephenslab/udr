# Compute the log-likelihood for the Ultimate Deconvolution model.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud <- function (X, w, U, V, version = c("Rcpp","R")) {
  m <- ncol(X)
  version <- match.arg(version)
  if (is.list(U))
    U <- simplify2array(U)
  if (is.matrix(V)) {
    if (version == "Rcpp")
      y <- loglik_ud_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_helper(X,w,U,V)
  } else {
    if (is.list(V))
      V <- simplify2array(V)
    if (version == "Rcpp")
      y <- loglik_ud2_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud2_helper(X,w,U,V)
  }
  return(y)
}

# This is the R implementation of the likelihood computation when the
# residual covariance V is the same for all samples.
loglik_ud_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (j in 1:k)
    y <- y + w[j] * dmvnorm(X,sigma = V + U[,,j])
  return(sum(log(y)))
}

# This is the R implementation of the likelihood computation when the
# residual covariance V is *not* the same for all samples.
loglik_ud2_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (i in 1:n)
    for (j in 1:k)
      y[i] <- y[i] + w[j] * dmvnorm(X[i,],sigma = V[,,i] + U[,,j])
  return(sum(log(y)))
}
