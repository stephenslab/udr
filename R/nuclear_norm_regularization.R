# This file implements nuclear norm penalty to covariance estimation.

#' Function to calculate the partial derivative w.r.t. one eigenvalue
#' of U. See equation (44) in https://www.overleaf.com/read/vrgwpskkhbpj.
grad_loglik_per_eigenval <- function(n, q, d, lambda, alpha){
  grad = n/(q+1)-n*d/(q+1)^2+lambda*alpha-lambda*(1-alpha)/q^2
  return(grad)
}

#' Function to regularize U by nuclear penalty from Chi and Lange (2014).
#' Here I assume V_j = I. U is the estimate from TED/ED/FA.
#' @param X: data matrix
#' @param lambda: the strength of penalty
#' @param alpha: control the trade-off between the two nuclear norm terms.
regularize_by_nuclear_penalty = function(X, lambda, alpha){
  n = nrow(X)
  m = ncol(X)
  S = t(X) %*% X/n
  evd = eigen(S)
  d = evd$values
  q = rep(0, m)
  for (i in 1:m){
    q[i] = uniroot(function(q) grad_loglik_per_eigenval(n, q, d[i],lambda,alpha),
                   c(0,1e6))$root
  }
  U = evd$vectors %*% diag(q) %*% t(evd$vectors)
  return(list(U = U, shrinked_eigenval = q))
}

#' Function to get alpha based on eq.(5) from Chi and Lange (2014).
#' This choice of alpha results in the prior mode being \hat\sigma*I.
get_alpha = function(X){
  m = ncol(X)
  n = nrow(X)
  S = t(X) %*% X/ n
  alpha = 1/(1+(sum(diag(S))/m)^2)
  return(alpha)
}
