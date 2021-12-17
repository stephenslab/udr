# Compute n x k log-likelihood matrix.
compute_loglik_matrix <- function (fit, version = c("Rcpp","R")) {
  version <- match.arg(version)
  
  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")
  
  # Get the prior (U) and residual (V) covariance matrices.
  U <- ulist2array(fit$U)
  V <- fit$V
  if (is.list(V))
    V <- list2array(V)

  # Compute the loglik matrix.
  if (is.matrix(V)) {

    # Perform the computations for the special case when the same
    # residual variance is used for all samples. Here the R
    # implementation is fast so no C++ implementation is provided.
    fit$loglik_matrix <- compute_loglik_matrix_iid(fit$X,U,V)
  } else {
      
    # Perform the computations for the more general case when the
    # residual variance is not necessarily the same for all samples.
    if (version == "R")
      fit$loglik_matrix <- compute_loglik_matrix_notiid(fit$X,U,V)
    else if (version == "Rcpp")
      fit$loglik_matrix <- compute_loglik_matrix_notiid_rcpp(fit$X,U,V)
  }

  # Add row and column names to the loglik matrix.
  rownames(fit$loglik_matrix) <- rownames(fit$X)
  colnames(fit$loglik_matrix) <- names(fit$U)

  # Output the updated fit.
  return(fit)
}

# This implements the calculations for compute_loglik_matrix for the
# special case when the same residual covariance matrix is used for
# all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_loglik_matrix_iid <- function (X, U, V) {
  n <- nrow(X)
  k <- dim(U)[3]
  loglik <- matrix(0,n,k)
  for (j in 1:k)
    loglik[,j] = dmvnorm(X,sigma = V + U[,,j],log = TRUE)
  return(loglik)
}

# This implements the calculations for compute_loglik_matrix for the
# more general case when the residual covariances are not necessaily
# the same for all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_loglik_matrix_notiid <- function (X, U, V) {
  n <- nrow(X)
  k <- dim(U)[3]
  loglik <- matrix(0,n,k)
  for (i in 1:n)
    for (j in 1:k)
      loglik[i,j] = dmvnorm(X[i,],sigma = V[,,i] + U[,,j],log = TRUE)
  return(loglik)
}
