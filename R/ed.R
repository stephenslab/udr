# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), for the special case when the
# samples are i.i.d. (same V).
update_prior_covariance_unconstrained_ed_iid <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed_iid(X,U$mat,p))

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed_iid.
update_prior_covariance_unconstrained_ed_iid_rcpp <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed_iid_rcpp(X,U$mat,p))

# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), allow for different residual
# variances V among the samples.
update_prior_covariance_unconstrained_ed_notiid <- function (X, U, V, p, minval) {
  if (is.matrix(V))
    mat <- update_prior_covariance_ed_iid(X,U$mat,V,p)
  else
    mat <- update_prior_covariance_ed_general(X,U$mat,V,p)
  return(update_prior_covariance_unconstrained_struct(U,mat))
}

# Update the prior covariance matrix (U) in the model x ~ N(0,U + I)
# using the update formula derived in Bovy et al (2011). Input p is a
# vector of weights associated with the rows of X.
ed_iid <- function (X, U, p) {
  m  <- ncol(X)
  I  <- diag(m)
  B  <- solve(U + I,U)
  X1 <- crossprod((sqrt(safenormalize(p))*X) %*% B)
  return(U - U %*% B + X1)
}

# Perform an M-step update for a prior covariance matrix (U) using the
# update formula derived in Bovy et al (2011), allowing for residual
# covariance matrices that differ among the data samples (rows of X).
# @param p is a vector of weights
# @param U is a matrix
# @param V is a 3-d array, in which V[,,j] is the covariance matrix
# for the jth observation
update_prior_covariance_ed_general <- function (X, U, V, p) {
  n  <- nrow(X)
  m  <- ncol(X)
  B  <- array(0,c(m,m,n))
  bb <- array(0,c(m,m,n))
  for (i in 1:n) {
    Tinv    <- solve(U + V[,,i])
    bb[,,i] <- tcrossprod(((sqrt(p[i])*U) %*% Tinv) %*% X[i,])
    B[,,i]  <- p[i]*(U - (U %*% Tinv) %*% U)
  }
  Unew <- (sliceSums(bb) + sliceSums(B))/sum(p)
  return(Unew)
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ed.
update_prior_covariance_rank1_ed_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_rank1_ed_rcpp is not yet implemented") 
}

