#' This file implements nuclear norm penalty to covariance estimation
#' under scaled model. x_j\sim N(0, \sigma^2\tilde U + I) and the penalty 
#' is the nuclear norm penalty function on \tilde U. 


#' Function to calculate the partial derivative w.r.t. one eigenvalue
#' of U. See equation (80) in https://www.overleaf.com/read/vrgwpskkhbpj.
#' @param val: the ith eigenvalue that we want to solve for
#' @param d: the ith eigenvalue of weighted empirical covariance matrix.
#' @param p: the weight vector for a component
#' @param lambda: the strength of penalty
#' @param alpha: control the trade-off between the two nuclear norm terms.
#' @param sigma2: the scalar on U.
grad_loglik_per_eigenval_scaled <- function(val, p, d, lambda, alpha, sigma2){
  grad = sum(p)/(val+1)-sum(p)*d/(val+1)^2+lambda*alpha/sigma2-lambda*(1-alpha)*sigma2/val^2
  return(grad)
}

#' Function to regularize U by nuclear penalty from Chi and Lange (2014).
#' Here I assume V_j = I.
#' @param X: data matrix of size $n$ by $R$.
#' @param S: weighted empirical covariance matrix
#' @param p: a vector of size $n$, containing the posterior weight for a certain 
#' component for each observation. 
#' @param lambda: the strength of penalty
#' @param alpha: a number that controls the trade-off between the two nuclear norm terms.
#' @param sigma2: the scalar on U.
regularize_by_nuclear_penalty_scaled = function(X, S, p, lambda, alpha, sigma2){
  n = nrow(X)
  m = ncol(X)
  evd = eigen(S)
  d = evd$values
  eigenval = rep(0, m)
  for (i in 1:m){
    eigenval[i] = uniroot(function(val) grad_loglik_per_eigenval_scaled(val, p, d[i],lambda, alpha, sigma2),
                   c(1e-6,1e6))$root
  }
  U = evd$vectors %*% diag(eigenval) %*% t(evd$vectors)
  return(list(U = U, shrinked_eigenval = eigenval))
}

#' Function to update sigma2 in M-step optimization given 
#' all the eigenvalues of U. See eq.(82) in the writeup. 
#' @param alpha: a number that controls the trade-off between the two nuclear norm terms.
#' @param eigenvals: the eigenvalues of U.
get_sigma2 = function(alpha, eigenvals){
  sigma2 = alpha/(1-alpha)*sum(eigenvals)/(sum(1/eigenvals))
  sigma2 = sqrt(sigma2)
  return(sigma2)
}

#' Function to update nuclear norm penalized covariance matrix.
#' @param X: data matrix of size $n$ by $R$.
#' @param w: a vector of length $K$ containing initial mixture component weights.
#' @param U: a 3d array contains k prior covariance matrix
#' @param V: residual covariance matrix, set to be identity matrix.
#' @param lambda: tuning paramter that controls the strength of penalty
em_fit_scaled_nuclear_norm_penalized_update <- function(X, w, U, V, lambda, maxiter, tol, tol.obj){
  # initialize progress to store progress at each iteration
  
  progress <- as.data.frame(matrix(0, maxiter,5))
  names(progress) <- c("iter","loglik", "log_posterior", "delta.w","delta.u")
  progress$iter <- 1:maxiter

  # Get the number of samples (n) and the number of mixture components (k)
  n <- nrow(X)
  k <- length(w)
  m <- ncol(X)
  alpha <- rep(0.5,k)
  s <- rep(1, k)
  log_posterior0 = -Inf
  Fc = c()
  for (iter in 1:maxiter){
    # param values at previous iteration
    w0  <- w
    U0  <- U
    # E-step: calculate posterior probabilities using the current mu and sigmas
    P <- compute_posterior_probs_iid(X, w, U, V)
    #f <- compute_F(X, w, U, V, P, lambda, alpha)
    #Fc = c(Fc, f)
    # M-step: 
    # update covariance matrix 
    for (j in 1:k){
      S = get_S(X, P[,j])
      res = regularize_by_nuclear_penalty_scaled(X, S, P[, j], lambda, alpha[j], s[j])
      U[,,j] = res$U
      s[j] = get_sigma2(alpha[j], res$shrinked_eigenval)
    }
    # update mixture weight
    w = colSums(P)/n
    #f <- compute_F(X, w, U, V, P, lambda, alpha)
    #Fc = c(Fc, f)
    # Compute log-posterior
    log_posterior = compute_log_posterior_nuclear_scaled(X, w, U, V, lambda, alpha, s)
    loglik = loglik_ud_iid_helper(X,w,U,V)
    dw <- max(abs(w - w0))
    dU <- max(abs(U - U0))
    dlogpost <- log_posterior - log_posterior0
    progress[iter,"loglik"]  <- loglik
    progress[iter, "log_posterior"] <- log_posterior
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    dparam  <- max(dw,dU)
    log_posterior0 <- log_posterior
    if (dparam < tol | dlogpost < log(1 + tol.obj))
      break
    }
    return(list(w = w, U = U, progress = progress[1:iter, ], sigma = s, Fc = Fc))
}


#' Function to compute F-function after each E-step/M-step.
#' See eq.(69) in write-up https://www.overleaf.com/read/vrgwpskkhbpj. 
#' @param X: data matrix of size $n$ by $R$.
#' @param w: a vector of length $K$ containing mixture component weights.
#' @param U: a 3d array contains k prior covariance matrix
#' @param V: residual covariance matrix, set to be identity matrix.
#' @param P: the weight matrix for each observation under each component.
#' @param lambda: tuning parameter that controls the strength of penalty
#' @param alpha: tuning parameter controls the strength of two penalty terms. 
compute_F = function(X, w, U, V, P, lambda, alpha){
  log_prior = 0
  weighted_pi = 0
  n = nrow(X)
  loglik_mat = matrix(0, nrow=K, ncol=n)
  for(k in 1:K){
    eigenval = eigen(U[,,k])$values
    loglik_mat[k,] <- t(dmvnorm(X,sigma = U[,,k] + V,log = TRUE))
    log_prior = log_prior -lambda/2*(alpha[k]*sum(eigenval) + (1-alpha[k])*sum(1/eigenval))
    weighted_pi = weighted_pi + sum(P[,k])*log(w[k])
  }
  f = sum(t(loglik_mat)*P) + log_prior + weighted_pi - sum(P*log(P))
  return(f)
}

#' Function to compute log-posterior. The log-posterior should increase each iteration.
#' See eq.(64) in write-up https://www.overleaf.com/read/vrgwpskkhbpj. 
#' @param X: data matrix of size $n$ by $R$.
#' @param w: a vector of length $K$ containing mixture component weights.
#' @param U: a 3d array contains k prior covariance matrix
#' @param V: residual covariance matrix, set to be identity matrix.
#' @param lambda: tuning parameter that controls the strength of penalty
#' @param alpha: tuning parameter controls the strength of two penalty terms. 
#' @param s: a list of scalars for each U.
compute_log_posterior_nuclear_scaled <- function(X, w, U, V, lambda, alpha, s){
  log_prior <- 0
  K <- length(w)
  for (k in 1:K){
    eigenval = eigen(U[,,k])$values
    log_prior = log_prior -lambda/2*(alpha[k]*sum(eigenval)/s[k] + (1-alpha[k])*s[k]*sum(1/eigenval))
  }
  loglik <- loglik_ud_iid_helper(X, w, U, V) 
  log_posterior <- loglik + log_prior
  return(log_posterior)
}
