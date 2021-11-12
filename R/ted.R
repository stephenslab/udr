# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al
# (2013). Note that input U is not used, and is included only for
# consistency with the other update functions.
update_prior_covariance_unconstrained_ted_iid <- function (X, U, V, p,
                                                           minval = 0) {
  mat <- ted(X,V,p,minval)
  return(update_prior_covariance_struct_unconstrained(U,mat))
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ted_idd.
update_prior_covariance_unconstrained_ted_iid_rcpp <- function (X, U, V, p,
                                                                minval = 0) {
  mat <- ted_rcpp(X,V,p,minval,r = ncol(X))
  return(update_prior_covariance_struct_unconstrained(U,mat))
}

# Perform an M-step update for a rank-1 prior covariance matrix U
# using the "eigenvalue truncation" technique. Note that input U is
# not used, and is included only for consistency with the other update
# functions.
update_prior_covariance_rank1_ted_iid <- function (X, U, V, p, minval = 0) {
  mat <- ted(X,V,p,minval = minval,r = 1)
  vec <- getrank1(mat)
  return(update_prior_covariance_struct_rank1(U,vec))
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ted_idd. (This is not yet implemented,
# so currently it just calls the R function.)
update_prior_covariance_rank1_ted_iid_rcpp <- function (X, U, V, p,
                                                        minval = 0) {
  mat <- ted_rcpp(X,V,p,minval = minval,r = 1)
  vec <- getrank1(mat)
  return(update_prior_covariance_struct_rank1(U,vec))
}

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ted_notiid <- function (X, U, V, p)
  stop("unconstrained.update = \"ted\" does not work for case when data ",
       "points are not i.i.d. (different Vs)")
update_prior_covariance_rank1_ted_notiid <- function (X,U,V,p)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_rank1_ted_notiid_rcpp <- function (X,U,V,p)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_unconstrained_ted_notiid <- function (X,U,V,p)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_unconstrained_ted_notiid_rcpp <- function (X,U,V,p)
  update_prior_covariance_ted_notiid(X,U,V,p)

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al (2013).
# Input p is a vector of weights associated with the rows of X. Input
# r specifies an optional constraint on U; when r < n, where U is an n
# x n matrix, at most r of the eigenvalues are positive in the updated
# matrix.
ted <- function (X, V, p, minval = 0, r = ncol(X)) {
    
  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  p <- safenormalize(p)
  R <- chol(V)
  T <- crossprod((sqrt(p)*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite, with at most r of its eigenvalues being
  # positive.
  U <- shrink_cov(T,minval,r)
  
  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}
