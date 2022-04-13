# This file implements nuclear norm penalty to covariance estimation.

#' Function to calculate the partial derivative w.r.t. one eigenvalue
#' of U. See equation (55) in https://www.overleaf.com/read/vrgwpskkhbpj.
grad_loglik_per_eigenval <- function(p, q, d, lambda, alpha){
  grad = sum(p)/(q+1)-sum(p)*d/(q+1)^2+lambda*alpha-lambda*(1-alpha)/q^2
  return(grad)
}

#' Function to regularize U by nuclear penalty from Chi and Lange (2014).
#' Here I assume V_j = I. U is the estimate from TED/ED/FA.
#' @param X: data matrix of size $n$ by $R$.
#' @param p: a vector of size $n$, containing the posterior weight for a certain 
#' component for each observation. 
#' @param lambda: the strength of penalty
#' @param alpha: control the trade-off between the two nuclear norm terms.
regularize_by_nuclear_penalty = function(X, p, lambda, alpha){
  n = nrow(X)
  m = ncol(X)
  S = get_S(X, p)
  evd = eigen(S)
  d = evd$values
  q = rep(0, m)
  for (i in 1:m){
    q[i] = uniroot(function(q) grad_loglik_per_eigenval(p, q, d[i],lambda,alpha),
                   c(0,1e6))$root
  }
  U = evd$vectors %*% diag(q) %*% t(evd$vectors)
  return(list(U = U, shrinked_eigenval = q))
}

#' Function to get alpha based on eq.(5) from Chi and Lange (2014).
#' This choice of alpha results in the prior mode being \hat\sigma*I.
#' @param X: data matrix of size $n$ by $R$.
#' @param p: a vector of size $n$, containing the posterior weight for a certain 
#' component for each observation. 
get_alpha = function(X, p){
  m = ncol(X)
  S = get_S(X, p)
  alpha = 1/(1+(sum(diag(S))/m)^2)
  return(alpha)
}


#' Function to calculate weighted empirical covariance matrix. Each data point 
#' is weighted by its posterior weight for a certain component.
#' @param X: data matrix of size $n$ by $R$.
#' @param p: a vector of size $n$, containing the posterior weight for a certain 
#' component for each observation. 
get_S = function(X, p){
  n = nrow(X)
  S = crossprod(sqrt(p)*X)/n # t(Xw) %*% Xw/n 
  return(S)
}