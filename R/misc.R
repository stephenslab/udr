## This function computes loglikelihood
loglik.compute <- function (X, w, T) {
  n <- nrow(X)
  k <- length(w)
  y <- rep(0,n)
  for (j in 1:k)
    y <- y + w[j] * dmvnorm(X, sigma = T[[j]])
  return(sum(log(y)))
}

## This function takes as input an array of unnormalized
## log-probabilities logw and returns normalized probabilities such
## that the sum is equal to 1.
normalizelogweights <- function (logw){
  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  c <- max(logw)
  w <- exp(logw - c)
  # Normalize the probabilities.
  return(w/sum(w))
}

## This is the function for "shrinking" the covariance matrix T to get
## $\hat T_k$. Setting eigenvalues <1 to 1+eps.  eps resolves the
## numerical issues, ensuring chol() works for the output matrices.
shrink.cov = function(T, eps){
  evd = eigen(T)
  shrink_eigen = ifelse(evd$values > 1, evd$values, 1+eps)
  T.new = tcrossprod(evd$vectors %*% diag(sqrt(shrink_eigen)))
  return(T.new)
}

# Returns TRUE if and only if all marginal covariances T = S + U are
# s.p.d. (symmetric positive definite).
verify.marginal.covariances <- function (U, S) {
  k   <- length(U)
  out <- TRUE  
  for (i in 1:k)
    out <- out & is.matrix(tryCatch(chol(S + U[[i]]),
                                    error = function (e) NULL))
  return(out)
}

# Randomly generate initial estimates of the prior covariance matrices
# U by computing the sample covariances of random subsets of the data.
#
#' @importFrom stats cov
generate.random.covariances <- function (X, k) {

  # Get the number of data samples (n) and their dimension (m).
  n <- nrow(X)
  m <- ncol(X)
  
  # Select the size of the random subsets.
  n0 <- max(20,m + 2)
  if (n0 >= n)
    stop("Cannot generate random initial estimates of prior covariance ",
         "matrices; more data points are needed")
  
  # Generate random covariance matrices by randomly selecting small
  # subsets of the data, and computing the sample covariance from
  # these subsets.
  U        <- vector("list",k)
  names(U) <- paste0("k",1:k)
  for (i in 1:k)
    U[[i]] <- cov(X[sample(n,n0),])
  return(U)
}
