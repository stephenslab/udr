#' @importFrom stats logLik
#'
#' @title Title Goes Here
#' 
#' @description Description goes here.
#'
#' @param object Describe input argument "fit" here.
#'
#' @param version Description of input argument "version" here.
#'
#' @param \dots Additional arguments (unused).
#' 
#' @return Describe the return value here.
#' 
#' @method logLik ud_fit
#' 
#' @export
#' 
logLik.ud_fit <- function (object, version = c("Rcpp","R"), ...) {
  if (!(is.list(object) & inherits(object,"ud_fit")))
    stop("Input argument \"object\" should be an object of class \"ud_fit\",",
         "such as the output of ud_init")
  version <- match.arg(version)
  return(loglik_ud(object$X,object$w,object$U,object$V,version))
}

# Compute the log-likelihood for the Ultimate Deconvolution model.
# Input argument U should either be a list of length k, in which each
# U[[i]]$mat is an m x m matrix, or an m x m x k array. Input argument
# V should either be an m x m matrix, a list of matrices of length n,
# or an m x m x n array, where n is the number of data points.
loglik_ud <- function (X, w, U, V, version = c("Rcpp","R")) {
  version <- match.arg(version)
  
  # Process input arguments U and V as needed.
  if (is.list(U))
    U <- ulist2array(U)
  if (is.list(V))
    V <- list2array(V)

  if (is.matrix(V)) {

    # Compute the likelihood in the case when the residual variance is
    # the same for all samples.
    if (version == "Rcpp")
      y <- loglik_ud_iid_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_iid_helper(X,w,U,V)
  } else {

    # Compute the log-likelihood in the case when the residual
    # variance is not necessarily the same for all samples.
    if (version == "Rcpp")
      y <- loglik_ud_general_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_general_helper(X,w,U,V)
  }
  return(y)
}

# Compute the log-likelihood when the residual covariance V is the
# same for all samples.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud_iid_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (j in 1:k)
    y <- y + w[j] * dmvnorm(X,sigma = V + U[,,j])
  return(sum(log(y)))
}

# Compute the log-likelihood when the residual covariance V is not
# necessarily the same for all samples.
#
#' @importFrom mvtnorm dmvnorm
loglik_ud_general_helper <- function (X, w, U, V) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (i in 1:n)
    for (j in 1:k)
      y[i] <- y[i] + w[j] * dmvnorm(X[i,],sigma = V[,,i] + U[,,j])
  return(sum(log(y)))
}
