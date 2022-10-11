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
    stop("Input argument \"object\" should be an object of class \"ud_fit\"")
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
      y <- loglik_ud_notiid_rcpp(X,w,U,V)
    else if (version == "R")
      y <- loglik_ud_notiid_helper(X,w,U,V)
  }
  return(y)
}

# Compute the log-likelihood when the residual covariance V is the
# same for all samples.
loglik_ud_iid_helper <- function (X, w, U, V) {
  K <- dim(U)[3] # number of mixture components 
  n <- nrow(X)
  loglik_mat = matrix(0, nrow=K, ncol=n)
  for(k in 1:K){
    loglik_mat[k,] <- t(dmvnorm(X,sigma = U[,,k] + V,log = TRUE))
  }
  return(sum(apply(loglik_mat+log(w),2,log_sum_exp)))
}


# Compute the log-likelihood when the residual covariance V is not
# necessarily the same for all samples.
loglik_ud_notiid_helper <- function (X, w, U, V) {
  n <- nrow(X)
  y <- rep(0,n)
  for (i in 1:n)
    y[i] <- ldmvnormmix(X[i,],w,U,V[,,i])
  return(sum(y))
}



# Function to compute log-penalty for one prior covariance in the transformed
# space if the transformation is applied. 
# @param U: prior covariance matrix (in the transformed space). 
# U has to be unconstrained covariance, otherwise makes no sense to penalize.
# @param sigma2: The scalar attached to U.
# @param n0: the penalty strength of inverse-Wishart prior used in ED updates.
# @param lambda: the penalty strength of nuclear norm function used in TED updates.
# @param S0: a positive definite matrix used as the parameter of inverse-Wishart distribution. 
# @param alpha: a tuning parameter used in nuclear norm penalty function. 
# Default of 0.5 is recommended. 
compute_penalty <- function(U, sigma2, n0 = 0, lambda = 0, S0, alpha = 0.5){
  penalty_iw <- 0
  penalty_nu <- 0
  if (n0 != 0)
    penalty_iw = ldiwishart(W = U/sigma2, n0, S0)
  if (lambda != 0){
    eigenval = eigen(U)$values
    penalty_nu = -lambda/2*(alpha*sum(eigenval)/sigma2 + (1-alpha)*sigma2*sum(1/eigenval))
  }
  log_penalty <- penalty_iw + penalty_nu
  return(log_penalty)
}

# Function to compute penalized log-likelihood.
# @param loglik: the log-likelihood calculated on fitted mixture model.
# @param logplt: log of the penality
compute_loglik_penalized <- function(loglik, logplt){
  loglik_penalized <- loglik + logplt
  return(loglik_penalized)
}


# Compute the log-density of the inverse Wishart at W with n - d - 1
# degrees of freedom and scale matrix n*S, ignoring terms that do not
# depend on X.
ldiwishart <- function (W, n, S)
  -n/2*(ldet(W) + tr(S %*% solve(W)))