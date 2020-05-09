# Compute the log-likelihood for the multivariate normal means model.
#
#' @importFrom mvtnorm dmvnorm
loglik_mvebnm <- function (X, w, U, S, version = c("Rcpp","R")) {
  m <- ncol(X)
  k <- length(U)
  version <- match.arg(version)
  if (version == "Rcpp") {
    U <- array(simplify2array(U),c(m,m,k))  
    y <- loglik_mvebnm_rcpp(X,w,U,S)
  } else if (version == "R")
    y <- loglik_mvebnm_helper(X,w,U,S)
  return(y)
}

# This is the R implementation of the likelihood computation.
loglik_mvebnm_helper <- function (X, w, U, S) {
  n <- nrow(X)
  k <- length(U)
  y <- rep(0,n)
  for (i in 1:k)
    y <- y + w[i] * dmvnorm(X,sigma = S + U[[i]])
  return(sum(log(y)))
}
