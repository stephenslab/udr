# Initialize the data structure for a scaled prior covariance matrix.
create_prior_covariance_struct_scaled <- function (X, U, s = 1) {
  rownames(U) <- colnames(X)
  colnames(U) <- colnames(X)
  Q <- get_mat_Q(U)
  U <- list(s = s, U0 = U, mat = U, Q = Q)
  attr(U,"covtype") <- "scaled"
  return(U)
}

# Initialize the data structure for a rank-1 prior covariance matrix.
create_prior_covariance_struct_rank1  <- function (X, U) {
  u <- getrank1(U)
  U <- list(vec = u,mat = tcrossprod(u))
  names(U$vec)    <- colnames(X)
  rownames(U$mat) <- colnames(X)
  colnames(U$mat) <- colnames(X)
  attr(U,"covtype") <- "rank1"
  return(U)
}

# Initialize the data structure for an unconstrained prior covariance
# matrix.
create_prior_covariance_struct_unconstrained <- function (X, U, s = 1) {
  rownames(U) <- colnames(X)
  colnames(U) <- colnames(X)
  U <- list(s = s, mat = U)
  attr(U,"covtype") <- "unconstrained"
  return(U)
}

# Update the data structure for a scaled prior covariance matrix.
# Input "U" is the current data structure, and "s" is the new estimate
# of the scaling factor. This function is used in the
# update_prior_covariance_scaled_* functions.
update_prior_covariance_struct_scaled <- function (U, s) {
  U$s   <- s
  U$mat <- s * U$U0
  return(U)
}

# Update the data structure for a rank-1 prior covariance matrix.
# Input "U" is the current data structure, and "vec" is a vector
# containing the new estimates, such that the new rank-1 matrix is
# tcrossprod(vec). This function is used in the
# update_prior_covariance_rank1_* functions.
update_prior_covariance_struct_rank1 <- function (U, vec) {
  mat           <- tcrossprod(vec)
  names(vec)    <- names(U$vec)
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$vec         <- vec
  U$mat         <- mat
  return(U)
}

# Update the data structure for an unconstrained prior covariance
# matrix. Input "U" is the current data structure, and "mat" is the
# newly estimated matrix. This function is used in the
# update_prior_covariance_unconstrained_* functions.
update_prior_covariance_struct_unconstrained <- function (U, mat, s) {
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$mat <- mat
  U$s <- s
  return(U)
}

# If Cov(x) = U, then cov(y) = t(A)*U*A, where y = t(A)*x. This
# function applies the affine transformation y = a*x to the prior
# covariance matrix U in which U is a scaled matrix.
transform_prior_covariance_struct_scaled <- function (U, A) {
  x <- rownames(U$U0)
  U$U0 <- t(A) %*% U$U0 %*% A
  U$mat <- U$s * U$U0
  U$Q <- get_mat_Q(U$U0)
  rownames(U$U0)  <- x
  colnames(U$U0)  <- x
  rownames(U$mat) <- x
  colnames(U$mat) <- x
  return(U)
}

# If Cov(x) = U, then cov(y) = t(A)*U*A, where y = t(A)*x. This
# function applies the affine transformation y = a*x to the prior
# covariance matrix U in which U is a rank-1 unconstrained matrix.
transform_prior_covariance_struct_rank1 <- function (U, A) {
  x <- names(U$vec)
  U$vec <- drop(crossprod(A,U$vec))
  U$mat <- tcrossprod(U$vec)
  names(U$vec)    <- x
  rownames(U$mat) <- x
  colnames(U$mat) <- x
  return(U)
}

# If Cov(x) = U, then cov(y) = t(A)*U*A, where y = t(A)*x. This
# function applies the affine transformation y = a*x to the prior
# covariance matrix U in which U is an unconstrained matrix.
transform_prior_covariance_struct_unconstrained <- function (U, A) {
  x <- rownames(U$mat)
  U$mat <- t(A) %*% U$mat %*% A
  rownames(U$mat) <- x
  colnames(U$mat) <- x
  return(U)
}
