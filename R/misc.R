# This function takes as input a matrix of unnormalized
# log-probabilities, and returns normalized probabilities such that
# each row sums to 1.
normalizelogweights <- function (W) {

  # Guard against underflow or overflow by adjusting the
  # log-probabilities so that the largest probability is 1.
  a <- apply(W,1,max)
  W <- exp(W - a)
  
  # Normalize the probabilities.
  return(W / rowSums(W))
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

# Returns TRUE if and only if the matrix is symmetric positive
# definite.
is.posdef <- function (X)
  is.matrix(tryCatch(chol(X),error = function (e) NULL))

# Returns TRUE if and only if the matrix is (symmetric) positive
# semi-definite.
is.semidef <- function (X, e = 1e-8) 
  all(eigen(X)$values > -e)

# Returns TRUE if and only if all prior covariances U are positive
# semi-definite, and all marginal covariances T = S + U are positive
# definite.
verify.prior.covariances <- function (U, S, e = 1e-8) {
  k   <- length(U)
  out <- TRUE  
  for (i in 1:k)
    out <- out & is.semidef(U[[i]]) & is.posdef(S + U[[i]])
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
  U <- vector("list",k)
  names(U) <- paste0("k",1:k)
  for (i in 1:k)
    U[[i]] <- cov(X[sample(n,n0),])
  return(U)
}
