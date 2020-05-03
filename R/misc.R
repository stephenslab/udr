## This function computes loglikelihood
loglik.compute <- function (X, w, T) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (j in 1:k)
    y <- y + w[j] * dmvnorm(X, sigma = T[[j]])
  return(sum(log(y)))
}



## This function takes as input an array of unnormalized log-probabilities logw
## and returns normalized probabilities such that the sum is equal to 1.
normalizelogweights <- function (logw){
  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  c <- max(logw)
  w <- exp(logw - c)
  # Normalize the probabilities.
  return(w/sum(w))
}


## This is the function for "shrinking" the covariance matrix T to get $\hat T_k$. Setting eigenvalues <1 to 1+eps.
## eps resolves the numerical issues, ensuring chol() works for the output matrices.
shrink.cov = function(T, eps){
  evd = eigen(T)
  shrink_eigen = ifelse(evd$values > 1, evd$values, 1+eps)
  T.new = tcrossprod(evd$vectors %*% diag(sqrt(shrink_eigen)))
  return(T.new)
}

