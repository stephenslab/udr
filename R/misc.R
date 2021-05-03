# Convert an n x m x k array into a list of n x m matrices.
array2list <- function (x) {
  k <- dim(x)[3]
  y <- vector("list",k)
  for (i in 1:k)
    y[[i]] <- x[,,i]
  return(y)
}

# Convert a list of n x m matrices to an n x m x k array.
list2array <- function (x) {
  n <- nrow(x[[1]])
  m <- ncol(x[[2]])
  k <- length(x)
  return(array(simplify2array(x),c(n,m,k)))
}

# Convert a list of m x m prior covariance matrices ("U" matrices) to
# an m x m x k array, In the input, the m x m matrices are stored in
# x[[i]]$mat for each list element i.
ulist2array <- function (x)
  list2array(lapply(x,function (e) "[["(e,"mat")))

# Returns TRUE if and only if the matrix is symmetric positive
# definite.
isposdef <- function (X)
  is.matrix(tryCatch(chol(X),error = function (e) NULL))

# Returns TRUE if and only if the matrix is (symmetric) positive
# semi-definite.
issemidef <- function (X, minval = -1e-15) 
  all(eigen(X)$values > minval)

# For symmetric semi-definite matrix X, return u such that
# tcrossprod(u) is the nearest rank-1 matrix.
getrank1 <- function (X) {
  out <- eigen(X)
  return(with(out,sqrt(values[1]) * vectors[,1]))
}

# Return a matrix containing the sums over the "slices".
sliceSums <- function (x)
  rowSums(x,dims = 2)

# Find the n x n matrix U + I that best approximates T satisfying the
# constraint that U is positive definite. This is achieved by setting
# any eigenvalues of T less than 1 to 1 + minval, or, equivalently,
# setting any eigenvalues of U less than zero to be minval. The output
# is a positive definite matrix, U. When r < n, U is additionally
# constrained so that at most r of its eigenvalues are positive.
shrink_cov <- function (T, minval = 0, r = nrow(T)) {
  n <- nrow(T)
  r <- min(r,n)
  out <- eigen(T)
  d <- out$values
  d <- pmax(d-1,minval)
  if (r < n)
    d[seq(r+1,n)] <- minval
  return(tcrossprod(out$vectors %*% diag(sqrt(d))))
}

# Output y = x/sum(x), but take care of the special case when all the
# entries are zero, in which case return the vector of all 1/n, where
# n = length(x).
safenormalize <- function (x) {
  n <- length(x)
  if (sum(x) <= 0)
    y <- rep(1/n,n)
  else
    y <- x/sum(x)
  return(y)
}

# Compute the softmax of each row of W in a way that guards against
# numerical underflow or overflow. The return value is a matrix of the
# same dimension in which the entries in each row sum to 1.
softmax <- function (W) {
  a <- apply(W,1,max)
  W <- exp(W - a)
  return(W/rowSums(W))
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
  drop(rWishart(1,max(4,m),diag(m)))
