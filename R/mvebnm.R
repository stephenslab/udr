# This file contains the code for “TEEM", truncated eigenvalue
# expectation-maximization algorithm.  Major modifications comparing
# to EM.R, 1) Input exactly matches ed.R; 2) eps(epsilon) added to
# resolve numerical issues in estimation.


## This is the main function for running the EM algorithm to fit the
## "extreme deconvolution" mixture model.

#' @title Title Goes Here.
#'
#' @description Description goes here.
#'
#' @param X Describe input argument "X" here.
#'
#' @param w Describe input argument "U" here.
#'
#' @param S Describe input argument "S" here.
#' 
#' @references
#'
#' Add references here.
#' 
#' @useDynLib mvebnm
#'
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
mvebnm <- function (X, w, U, S, eps, maxiter, tol = 1e-6, verbose = TRUE) {
    
  # Input parameters:
  # X: n by m data matrix, n is number of samples, m is the number of outcomes per sample.
  # w: input weights for k mixture components, k by 1 vector.
  # U: initial estimates for covariance matrices, list of m by m matrices.
  # S: T = U + S, where T is the covariance matrix for zscore_hat. And in TEEM, S = I(identity matrix)
  # eps: we usually specify a small number for eps, say 1e-8, used to resolve numerical issues.
  # maxiter: maximum number of iterations
  # tol: criteria for convergence

  # initialize progress to store progress at each iteration
  progress = data.frame(iter = 1:maxiter,
                        obj  = rep(0,maxiter),
                        maxd = rep(0,maxiter))
                        
  n <- nrow(X)
  k <- length(w)
  m <- ncol(X)
  T = c()
  
  # get a list of T. (T = U + S)
  for (i in 1:k)
    T[[i]] = U[[i]] + S

  if (verbose)
    cat("iter         objective max.diff\n")

  for (iter in 1:maxiter){
      
    # store parameters and likelihood in the previous step
    w0  <- w
    T0  <- T

    # E-step: calculate posterior probabilities using the current mu and sigmas
    logP = matrix(0, nrow = n, ncol = k)
    for (j in 1:k){
      # log-density
      logP[,j] = log(w[j]) + dmvnorm(X, sigma = T[[j]], log = TRUE)
    }
    P = t(apply(logP, 1, function(x) normalizelogweights(x)))

    # M-step:
    # update covariance matrix with constraints
    for (j in 1:k){
      T[[j]] = t(X)%*%(P[,j]*(X))/sum(P[,j])
      T[[j]] = shrink.cov(T[[j]], eps)
    }

    # update mixture weight
    w = colSums(P)/n

    # Compute log-likelihood at the current estimates
    f = loglik.compute(X, w, T)
    
    d = max(abs(w - w0))
    progress[iter,"obj"]  <- f
    progress[iter,"maxd"] <- d


    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e\n",iter,f,d))

    if (d < tol)
      break
  }
  
  for (i in 1:k){
    U[[i]] = T[[i]] - S
  }
  return(list(w = w, U = U, progress = progress[1:iter, ]))
}
