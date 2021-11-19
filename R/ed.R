# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), for the special case when V =
# I for all samples; that is, the model is x ~ N(0,U + I).
update_prior_covariance_unconstrained_ed_iid <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed_iid(X,U$mat,p))

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed_iid.
update_prior_covariance_unconstrained_ed_iid_rcpp <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed_iid_rcpp(X,U$mat,p))

# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), allow for different residual
# variances V among the samples.
update_prior_covariance_unconstrained_ed_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed(X,U$mat,V,p))

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ed_notiid. (The C++ version has not yet been
# implemented, so for now we simply all the R implementation.)
update_prior_covariance_unconstrained_ed_notiid_rcpp <- function (X, U, V, p,
                                                                  ...)
  update_prior_covariance_struct_unconstrained(U,ed(X,U$mat,V,p))

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ed_invalid <- function (X, U, V, p, ...)
  stop("control$scale.update = \"ed\" and control$rank1.update = \"ed\" ",
       "are not valid choices")
update_prior_covariance_scaled_ed_iid         <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_iid_rcpp    <- function (X, U, p,...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_iid          <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_iid_rcpp     <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_notiid      <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_notiid_rcpp <- function (X, U, V, p,...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_notiid       <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_notiid_rcpp  <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()

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
ed <- function (X, U, V, p) {
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


