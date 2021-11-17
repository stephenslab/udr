# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al
# (2013). Note that input U is not used, and is included only for
# consistency with the other update functions.
update_prior_covariance_unconstrained_ted_iid <-
  function (Y, U, R, p, minval = 0, ...) {
  mat <- ted(Y,p,minval)
  mat <- t(R) %*% mat %*% R
  return(update_prior_covariance_struct_unconstrained(U,mat))
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ted_iid.
update_prior_covariance_unconstrained_ted_iid_rcpp <-
  function (Y, U, R, p, minval = 0, ...) {
  mat <- ted_rcpp(Y,p,minval,r = ncol(Y))
  mat <- t(R) %*% mat %*% R
  return(update_prior_covariance_struct_unconstrained(U,mat))
}

# Perform an M-step update for a rank-1 prior covariance matrix U
# using the "eigenvalue truncation" technique. Note that input U is
# not used, and is included only for consistency with the other update
# functions.
update_prior_covariance_rank1_ted_iid <-
  function (Y, U, R, p, minval = 0, ...) {
  mat <- ted(Y,p,minval = minval,r = 1)
  mat <- t(R) %*% mat %*% R
  return(update_prior_covariance_struct_rank1(U,getrank1(mat)))
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ted_iid.
update_prior_covariance_rank1_ted_iid_rcpp <-
  function (Y, U, R, p, minval = 0, ...) {
  mat <- ted_rcpp(Y,p,minval = minval,r = 1)
  mat <- t(R) %*% mat %*% R
  return(update_prior_covariance_struct_rank1(U,getrank1(mat)))
}

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ted_notiid <- function (X, U, V, p, ...)
  stop("unconstrained.update = \"ted\" does not work for case when data ",
       "points are not i.i.d. (different Vs)")
update_prior_covariance_rank1_ted_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(Y,U,R,p)
update_prior_covariance_rank1_ted_notiid_rcpp <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(Y,U,R,p)
update_prior_covariance_unconstrained_ted_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(Y,U,R,p)
update_prior_covariance_unconstrained_ted_notiid_rcpp <- function (X,U,V,p,...)
  update_prior_covariance_ted_notiid(Y,U,R,p)

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al
# (2013). (This update can also be viewed as computing
# maximum-likelihood estimates of the parameters in the probabilistic
# PCA model.) Input p is a vector of weights associated with the rows
# of Y. Input r specifies an optional constraint on U; when r < n,
# where U is an n x n matrix, at most r of the eigenvalues are
# positive in the updated matrix.
ted <- function (Y, p, minval = 0, r = ncol(Y)) {
    
  # Account for the weights by scaling the rows of Y.
  T <- crossprod(sqrt(safenormalize(p))*Y)
  
  # Find U maximizing the likelihood under the model y ~ N(0,I + U),
  # subject to U being a symmetric positive definite (s.p.d.) matrix,
  # and with the additional constraint that at most r of its
  # eigenvalues are positive.
  return(shrink_cov(T,minval,r))
}
