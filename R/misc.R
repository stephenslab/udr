# Convert an n x m x k array into a list of n x m matrices.
array2list <- function (x) {
  k <- dim(x)[3]
  y <- vector("list",k)
  for (i in 1:k)
    y[[i]] <- x[,,i]
  return(y)
}

# Returns TRUE if and only if the matrix is symmetric positive
# definite.
isposdef <- function (X)
  is.matrix(tryCatch(chol(X),error = function (e) NULL))

# Returns TRUE if and only if the matrix is (symmetric) positive
# semi-definite.
issemidef <- function (X, e = 1e-15) 
  all(eigen(X)$values > -e)

# For symmetric semi-definite matrix X, get the nearest rank-1 matrix.
getrank1 <- function (X) {
  out <- eigen(X)
  out$values[-1] <- 0
  return(with(out,vectors %*% diag(values) %*% t(vectors)))
}
 
# Compute the softmax of each row of W in a way that guards against
# numerical underflow or overflow. The return value is a matrix of the
# same dimension in which the entries in each row sum to 1.
softmax <- function (W) {
  a <- apply(W,1,max)
  W <- exp(W - a)
  return(W / rowSums(W))
}

# "Shrink" matrix T = U + I; that is, find the "best" matrix T
# satisfying the constraint that U is positive definite. This is
# achieved by setting any eigenvalues of T less than 1 to 1 + minval,
# or, equivalently, setting any eigenvalues of U less than 0 to be
# minval.  The output is a positive definite matrix, U.
shrink_cov <- function (T, minval) {
  if (length(T) == 1)

    # Handle univariate case (m = 1).
    U <- pmax(T - 1,minval)  
  else {
    out <- eigen(T)
    d   <- out$values
    d   <- pmax(d - 1,minval)
    U   <- tcrossprod(out$vectors %*% diag(sqrt(d)))
  }
  return(U)
}

# Returns TRUE if and only if all prior covariances U are positive
# semi-definite, and all marginal covariances T = V + U are positive
# definite.
verify_prior_covariances <- function (U, V, e = 1e-15) {
  k   <- length(U)
  out <- TRUE  
  for (i in 1:k)
    out <- out & issemidef(U[[i]]) & isposdef(V + U[[i]])
  return(out)
}

# Randomly generate an m x m symmetric rank-1 matrix.
#
#' @importFrom stats rnorm
sim_rank1 <- function (m)
  tcrossprod(rnorm(m))

# Randomly generate an m x m (unconstrained) symmetric matrix.
#
#' @importFrom stats rWishart
sim_unconstrained <- function (m)
  drop(rWishart(1,4,diag(m)))
