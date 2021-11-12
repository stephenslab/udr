# Update the scaling factor for a prior canonical covariance (U) matrix.
# @param U0 A (fixed) covariance matrix.
# @param V A covariance matrix.
# @return An integer scalar
#
#' @importFrom stats uniroot
update_prior_covariance_scaled_iid <- function (X, U0, V, p, minval) {

  # Transform data using the trick
  # Uhat = R^{-T}*U*R^{-1}
  # Xhat = X*R^{-1}
  R    <- chol(V)
  Uhat <- solve(t(R),U0) %*% solve(R)
  Xhat <- X %*% solve(R)
    
  # Eigenvalue decomposition based on transformed U.
  evd <- eigen(Uhat)
  lambdas <- ifelse(evd$values < minval,minval,evd$values)

  Y <- t(Xhat %*% evd$vectors)
  dev_zero <- grad_loglik_scale_factor(0, p, Y, lambdas) 
  if (dev_zero <= 0){
    return(0)
  }else{
    return(uniroot(function (s) grad_loglik_scale_factor(s,p,Y,lambdas),
                   c(0,1e6))$root)
  }
}

# Function for 1-d search of s value based on eq. (20) in the write-up.
# @param p Vector of weights.
# @param s The scaling factor for one component we aim to search for
# @param Y The transformed data
# @param lambdas Eigenvalues of U.
# @return A function of the scalar s.
grad_loglik_scale_factor <- function (s, p, Y, lambdas)
  sum(p*apply(Y,2,function(y) sum(lambdas*y^2/((s*lambdas + 1)^2)) -
                sum(lambdas/(s*lambdas + 1))))
