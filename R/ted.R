# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al
# (2013). Note that input U is not used, and is included only for
# consistency with the other update functions.
update_prior_covariance_unconstrained_ted_iid <- function (X, U, p,
                                                           minval = 0, lambda, penalty.type, ...){
  if (lambda == 0){
    return(update_prior_covariance_struct_unconstrained(U, ted(X,p,minval), s = 1))
  }
  else{
    S <- get_S(X, p)
    # Update U and store penalized eigenvalues
    if (penalty.type == "nu"){
      res <- ted_penalized_nu(X, S, p, lambda, alpha = 0.5, U$s) 
      sigma2 <- get_sigma2_nu(alpha = 0.5, res$shrinked_eigenval)
    }
    if (penalty.type == "iw"){
      res <- ted_penalized_iw(X, S, p, lambda, U$s)
      sigma2 <- get_sigma2_iw(res$shrinked_eigenval)
    }
  return(update_prior_covariance_struct_unconstrained(U, res$U, sigma2))
  }
}
  

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ted_iid.
update_prior_covariance_unconstrained_ted_iid_rcpp <-
  function (X, U, p, minval = 0, ...)
  update_prior_covariance_struct_unconstrained(U,ted_rcpp(X,p,minval,
                                                          r = ncol(X)))

# Perform an M-step update for a rank-1 prior covariance matrix U
# using the "eigenvalue truncation" technique. Note that input U is
# not used, and is included only for consistency with the other update
# functions.
update_prior_covariance_rank1_ted_iid <- function (X, U, p, minval = 0, ...)
  update_prior_covariance_struct_rank1(U,getrank1(ted(X,p,minval = minval,
                                                      r = 1)))

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ted_iid.
update_prior_covariance_rank1_ted_iid_rcpp <- function (X, U, p, minval = 0,
                                                        ...)
  update_prior_covariance_struct_rank1(U,getrank1(ted_rcpp(X,p,minval = minval,
                                                           r = 1)))

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ted_invalid <- function()
  stop("control$unconstrained.update = \"ted\" does not work for case when ",
       "data points are not i.i.d. (different Vs)")
update_prior_covariance_rank1_ted_notiid              <- function (X,U,V,p,...)
  update_prior_covariance_ted_invalid()
update_prior_covariance_rank1_ted_notiid_rcpp         <- function (X,U,V,p,...)
  update_prior_covariance_ted_invalid()
update_prior_covariance_unconstrained_ted_notiid      <- function (X,U,V,p,...)
  update_prior_covariance_ted_invalid()
update_prior_covariance_unconstrained_ted_notiid_rcpp <- function (X,U,V,p,...)
  update_prior_covariance_ted_invalid()

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al (2013).
# (This update can also be viewed as computing maximum-likelihood
# estimates of the parameters in the probabilistic PCA model.) Input p
# is a vector of weights associated with the rows of X. Input r
# specifies an optional constraint on U; when r < n, where U is an n x
# n matrix, at most r of the eigenvalues are positive in the updated
# matrix.
ted <- function (X, p, minval = 0, r = ncol(X)) {
    
  # Account for the weights by scaling the rows of X.
  T <- crossprod(sqrt(safenormalize(p))*X)
  
  # Find U maximizing the likelihood under the model x ~ N(0,U + I),
  # subject to U being a symmetric positive definite (s.p.d.) matrix,
  # and with the additional constraint that at most r of its
  # eigenvalues are positive.
  return(shrink_cov(T,minval,r))
}




# The following functions implement nuclear norm penalty to covariance
# estimation under scaled model, x_j ~ N(0, \sigma^2\tilde U + I).
# The nuclear norm penalty is specified on \tilde{U}.

# Function to compute the negative of objective function to be minimized. 
# @param val: the ith eigenvalue that we want to solve for
# @param d: the ith eigenvalue of weighted empirical covariance matrix.
# @param p: the weight vector for a component
# @param lambda: the strength of penalty
# @param alpha: control the trade-off between the two nuclear norm terms.
# @param sigma2: the scalar on U.
loglik_neg_per_eigenval_nu <- function(val, p, d, lambda, alpha, sigma2){
  obj = sum(p)*log(val+1) + sum(p)*d/(val + 1) + lambda * (alpha/sigma2*val+(1-alpha)*sigma2/val)
  return(obj/2)
}

# Function to calculate the partial derivative w.r.t. one eigenvalue
# of U. See equation (80) in https://www.overleaf.com/read/vrgwpskkhbpj.
# @param val: the ith eigenvalue that we want to solve for
# @param d: the ith eigenvalue of weighted empirical covariance matrix.
# @param p: the weight vector for a component
# @param lambda: the strength of penalty
# @param alpha: control the trade-off between the two nuclear norm terms.
# @param sigma2: the scalar on U.
grad_loglik_per_eigenval_nu <- function(val, p, d, lambda, alpha, sigma2){
  grad = sum(p)/(val+1)-sum(p)*d/(val+1)^2+lambda*alpha/sigma2-lambda*(1-alpha)*sigma2/val^2
  return(grad)
}


# Function to regularize U by nuclear penalty from Chi and Lange (2014).
# Here I assume V_j = I.
# @param X: data matrix of size $n$ by $R$.
# @param S: weighted empirical covariance matrix
# @param p: a vector of size $n$, containing the posterior weight for a certain 
# component for each observation. 
# @param lambda: the strength of penalty
# @param alpha: a number that controls the trade-off between the two nuclear norm terms.
# @param sigma2: the scalar on U.
ted_penalized_nu <- function(X, S, p, lambda, alpha, sigma2){
  n = nrow(X)
  m = ncol(X)
  evd = eigen(S)
  d = evd$values
  eigenval = rep(0, m)
  for (i in 1:m){
    eigenval[i] = optim(par = 1, fn = loglik_neg_per_eigenval_nu, gr = grad_loglik_per_eigenval_nu,
                        p, d[i], lambda, alpha, sigma2, method = "L-BFGS-B", lower = 1e-8, upper = 1e6)$par
  }
  U = evd$vectors %*% diag(eigenval) %*% t(evd$vectors)
  return(list(U = U, shrinked_eigenval = eigenval))
}

# Function to update sigma2 in M-step optimization given 
# all the eigenvalues of U. See eq.(82) in the writeup. 
# @param alpha: a number that controls the trade-off between the two nuclear norm terms.
# @param eigenvals: the eigenvalues of U.
get_sigma2_nu = function(alpha, eigenvals){
  sigma2 = alpha/(1-alpha)*sum(eigenvals)/(sum(1/eigenvals))
  sigma2 = sqrt(sigma2)
  return(sigma2)
}

# Function to calculate weighted empirical covariance matrix. Each data point 
# is weighted by its posterior weight for a certain component.
# @param X: data matrix of size $n$ by $R$.
# @param p: a vector of size $n$, containing the posterior weight for a certain 
# component for each observation. 
get_S = function(X, p){
  S = crossprod(sqrt(p)*X)/sum(p) 
  return(S)
}


# Function to calculate the partial derivative in IW penalty, w.r.t. one eigenvalue
# of U. See (8.70) in udr paper: https://www.overleaf.com/read/vdnzcbpgmgck
# @param val: the ith eigenvalue that we want to solve for
# @param d: the ith eigenvalue of weighted empirical covariance matrix.
# @param p: the weight vector for a component
# @param lambda: the strength of penalty
# @param sigma2: the scalar on U.
grad_loglik_per_eigenval_scaled_iw <- function(val, p, d, lambda, sigma2){
  grad = sum(p)/(val+1)-sum(p)*d/(val+1)^2+lambda*(1/val-sigma2/(val^2))
  return(grad)
}

# IW penalty: the negative of objective function to be minimized. 
# @param val: the ith eigenvalue that we want to solve for
# @param d: the ith eigenvalue of weighted empirical covariance matrix.
# @param p: the weight vector for a component
# @param lambda: the strength of penalty
# @param sigma2: the scalar on U.
loglik_neg_per_eigenval_scaled_iw <- function(val, p, d, lambda, sigma2){
  obj = sum(p)*log(val+1) + sum(p)*d/(val + 1) + lambda * (log(val)+sigma2/val)
  return(obj/2)
}

# IW penalty: function to update sigma2 in M-step optimization given 
# all the eigenvalues of U. See eq.(8.72) in the udr paper Appendix. 
# @param eigenvals: the eigenvalues of U.
get_sigma2_iw = function(eigenvals){
  R = length(eigenvals)
  sigma2 = R/sum(1/eigenvals)
  return(sigma2)
}

#' Function to regularize U by IW penalty. Here I assume V_j = I.
# @param X: data matrix of size $n$ by $R$.
# @param S: weighted empirical covariance matrix
# @param p: a vector of size $n$, containing the posterior weight for a certain 
# component for each observation. 
# @param lambda: the strength of penalty
# @param sigma2: the scalar on U.
ted_penalized_iw <- function(X, S, p, lambda, sigma2){
  n = nrow(X)
  m = ncol(X)
  evd = eigen(S)
  d = evd$values
  eigenval = rep(0, m)
  for (i in 1:m){
    eigenval[i] = optim(par = 1, fn = loglik_neg_per_eigenval_scaled_iw, gr = grad_loglik_per_eigenval_scaled_iw,
                        p, d[i], lambda, sigma2, method = "L-BFGS-B", lower = 1e-8, upper = 1e6)$par
  }
  U = evd$vectors %*% diag(eigenval) %*% t(evd$vectors)
  return(list(U = U, shrinked_eigenval = eigenval))
}

