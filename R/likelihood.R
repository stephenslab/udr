#' @title Ultimate Deconvolution Model Likelihoods
#' 
#' @description Compute the log-likelihood for the Ultimate
#' Deconvolution model.
#'
#' @param object An Ultimate Deconvolution model fit. Typically,
#'   this will be an output of \code{\link{ud_init}} or \code{ud_fit}.
#'
#' @param version When \code{version == "R"}, the computations are
#'   performed entirely in R; when \code{version == "Rcpp"}, an Rcpp
#'   implementation is used.
#'
#' @param \dots Additional arguments (unused).
#' 
#' @return A number giving the log-likelihood for the model.
#'
#' @seealso \code{\link{ud_init}}, \code{\link{ud_fit}}
#' 
#' @method logLik ud_fit
#' 
#' @importFrom stats logLik
#'
#' @export
#' 
logLik.ud_fit <- function (object, version = c("Rcpp","R"), ...) {
  version <- match.arg(version)
  if (!(is.list(object) & inherits(object,"ud_fit")))
    stop("Input argument \"object\" should be an object of class \"ud_fit\",",
         "such as the output of ud_init")
  out <- loglik_ud(object$X,object$w,object$U,object$V,version)
  class(out) <- "logLik"
  attr(out,"df") <- as.numeric(NA)
  return(out)
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
