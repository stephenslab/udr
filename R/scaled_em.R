
# Perform an M-step update for a scaled prior covariance (U) by directly
# maximizing the log-likelihood, for the special case when V = I for all samples; 
# that is, the model is x ~ N(0,sU + I).
update_prior_covariance_scaled_em_iid <- function (X, U, p, ...){
  update_prior_covariance_struct_scaled(U, em_scaled_iid(X, U$U0, p, minval = 0))
}


# Perform an M-step update for a scaled prior covariance (U) by directly
# maximizing the log-likelihood, allowing for different 
# residual variances V among the samples.
update_prior_covariance_scaled_em_notiid <- function (X, U, V, p, ...){
  update_prior_covariance_struct_scaled(U, em_scaled_notiid(X, U$U0, V, p, U$s))
}


# Update the scaling factor for a prior canonical covariance (U) matrix in iid case. 
#' @param U0 A (fixed) covariance matrix.
#' @param V A covariance matrix.
#' @return A scalar
#' @importFrom stats uniroot
em_scaled_iid <- function(X, U, p, minval){
  evd <- eigen(U)
  lambdas <- ifelse(evd$values < minval,minval,evd$values)
  Y <- t(X %*% evd$vectors)
  dev_zero <- grad_loglik_scaled_iid(0, p, Y, lambdas) 
  if (dev_zero <= 0){
    return(0)
  }else{
    return(uniroot(function (s) grad_loglik_scaled_iid(s,p,Y,lambdas),
                   c(0,1e6))$root)
  }
}


# Update the scaling factor for a prior canonical covariance (U) matrix in non-iid case. 
em_scaled_notiid <- function(X, U, V, p, s){
  s <- optim(par = 10, fn = weighted_loglik_negative, gr = grad_loglik_scaled_notiid,
            X = X, U = U, V = V, p = p , method = "L-BFGS-B", lower = 0, upper = 1e6)$par
  return(s)
}


# Function to calculate the gradient in scaled-EM model in iid case. 
#' @param s The scaling factor we aim to search for
#' @param p Vector of weights.
#' @param Y The transformed data
#' @param lambdas Eigenvalues of U.
#' @return A function of the scalar s.
grad_loglik_scaled_iid <- function (s, p, Y, lambdas)
  sum(p*apply(Y,2,function(y) sum(lambdas*y^2/((s*lambdas + 1)^2)) -
                sum(lambdas/(s*lambdas + 1))))


# Function to calculate the negative weighted loglikelihood 
# for optim() in non-iid case. 
#' @return A function of the scalar s.
weighted_loglik_negative <- function(s, X, U, V, p){
  weighted_logliks <- 0
  n <- nrow(X)
  for (i in 1:n){
    weighted_logliks <-  weighted_logliks + p[i]*dmvnorm(X[i, ], sigma = s*U+V[,,i], log = TRUE)
  }
  return(-weighted_logliks)
}



# Function to calculate the gradient of negative weighted log-likelihood
# in scaled-EM model in non-iid case for optim(). 
#' @return A function of the scalar s.
grad_loglik_scaled_notiid <- function (s, X, U, V, p){ 
  m <- ncol(U)
  n <- nrow(X)
  B = array(0, dim = c(m, m, n))
  Y = matrix(0, ncol = m, nrow = n)
  
  for (i in 1:n){
    B[,,i] = solve(s*U+V[,,i])
    Y[i, ] <- t(X[i, ]) %*% B[,,i] 
  }
  Yw <- sqrt(p) * Y
  trsum <- sum(diag(sliceSums(B*p) %*% U))
  gradient <- -trsum + sum(diag(Yw %*% U %*% t(Yw)))
  return(-gradient)
}









