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

# "Shrink" matrix T = U + I; that is, find the "best" matrix T
# satisfying the constraint that U is positive definite. This is
# achieved by setting any eigenvalues of T less than 1 to 1+e or,
# equivalently, setting any eigenvalues of U less than 0 to be e.
# The return value is the positive definite matrix U.
shrink.cov <- function (T, e) {
  out    <- eigen(T)
  values <- out$values
  values <- pmax(values - 1,e)
  U      <- tcrossprod(out$vectors %*% diag(sqrt(values)))
  return(U)
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
