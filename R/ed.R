# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011).
update_prior_covariance_unconstrained_ed_iid <- function (X, U, V, p, minval) {
  if (is.matrix(V))
    mat <- update_prior_covariance_ed_iid(X,U$mat,V,p)
  else
    mat <- update_prior_covariance_ed_general(X,U$mat,V,p)
  return(update_prior_covariance_unconstrained_struct(U,mat))
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed.
update_prior_covariance_unconstrained_ed_rcpp <- function (X, U, V, p,
                                                           minval) {
  if (is.matrix(V)) 
    mat <- update_prior_covariance_ed_iid_rcpp(X,U$mat,V,p)
  else
    stop("update_prior_covariance_unconstrained_ed_rcpp is not yet ",
         "implemented for case when data points are not i.i.d. (different Vs)")
  return(update_prior_covariance_unconstrained_struct(U,mat))
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

# Perform an M-step update for a prior covariance matrix (U) using the
# update formula derived in Bovy et al (2011), for the special case
# when the residual covariances V are the same for all data points.
# Input p is a vector of weights associated with the rows of X.
update_prior_covariance_ed_iid <- function (X, U, V, p) {
  p <- safenormalize(p)
  T <- U + V
  B <- solve(T,U)
  X1 <- crossprod((sqrt(p)*X) %*% B)
  return(U - U %*% B + X1)
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ed.
update_prior_covariance_rank1_ed_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_rank1_ed_rcpp is not yet implemented") 
}

