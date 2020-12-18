# Compute the log-likelihood for the Ultimate Deconvolution model.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud <- function (X, w, U, V, version = c("Rcpp","R")) {
  m <- ncol(X)
  k <- length(U)
  version <- match.arg(version)
  if (version == "Rcpp")
    y <- loglik_ud_rcpp(X,w,U,V)
  else if (version == "R")
    y <- loglik_ud_helper(X,w,U,V)
  return(y)
}

# This is the R implementation of the likelihood computation.
loglik_ud_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (i in 1:k)
    y <- y + w[i] * dmvnorm(X,sigma = V + U[,,i])
  return(sum(log(y)))
}
