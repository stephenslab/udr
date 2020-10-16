# Compute the log-likelihood for the multivariate normal means model.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud <- function (X, w, U, S, version = c("Rcpp","R")) {
  m <- ncol(X)
  k <- length(U)
  version <- match.arg(version)
  if (version == "Rcpp")
    y <- loglik_ud_rcpp(X,w,U,S)
  else if (version == "R")
    y <- loglik_ud_helper(X,w,U,S)
  return(y)
}

# This is the R implementation of the likelihood computation.
loglik_ud_helper <- function (X, w, U, S) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (i in 1:k)
    y <- y + w[i] * dmvnorm(X,sigma = S + U[,,i])
  return(sum(log(y)))
}
