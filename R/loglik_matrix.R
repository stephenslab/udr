#' @export
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
    # residual variance is used for all samples.
     if (version == "R")
      fit$loglik_matrix <- compute_loglik_matrix_iid(fit$X,U,V)
    else if (version == "Rcpp")
      fit$loglik_matrix <- compute_loglik_matrix_iid_rcpp(fit$X,U,V)
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
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  K <- dim(U)[3]

  # Compute the log-probabilities, stored in an n x k matrix.
  loglik_matrix <- matrix(0,n,K)
  for (k in 1:K)
    loglik_matrix[,k] = dmvnorm(X,sigma = V + U[,,k],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(loglik_matrix)
}

# This implements the calculations for compute_loglik_matrix for the
# more general case when the residual covariances are not necessaily
# the same for all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_loglik_matrix_notiid <- function (X, U, V) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  K <- dim(U)[3]

  # Compute the log-probabilities, stored in an n x k matrix.
  loglik_matrix <- matrix(0,n,K)
  for (i in 1:n)
    for (k in 1:K)
      loglik_matrix[i,k] = dmvnorm(X[i,],sigma = V[,,i] + U[,,k],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(loglik_matrix)
}


