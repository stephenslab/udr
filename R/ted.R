# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al
# (2013). Note that input U is not used, and is included only for
# consistency with the other update functions.
update_prior_covariance_unconstrained_ted_iid <- function (X, U, p,
                                                           minval = 0, ...)
  update_prior_covariance_struct_unconstrained(U,ted(X,p,minval))

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ted_iid.
update_prior_covariance_unconstrained_ted_iid_rcpp <-
  function (X, U, p, minval = 0, ...)
  update_prior_covariance_struct_unconstrained(U,ted_rcpp(X,p,minval,
                                                          r = ncol(X)))

# Perform an M-step update for a rank-1 prior covariance matrix U
# using the "eigenvalue truncation" technique. Note that input U is
# not used, and is included only for consistency with the other update
# functions.
update_prior_covariance_rank1_ted_iid <- function (X, U, p, minval = 0, ...)
  update_prior_covariance_struct_rank1(U,getrank1(ted(X,p,minval = minval,
                                                      r = 1)))

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ted_iid.
update_prior_covariance_rank1_ted_iid_rcpp <- function (X, U, p, minval = 0,
                                                        ...)
  update_prior_covariance_struct_rank1(U,getrank1(ted_rcpp(X,p,minval = minval,
                                                           r = 1)))

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ted_notiid <- function (X, U, V, p, ...)
  stop("unconstrained.update = \"ted\" does not work for case when data ",
       "points are not i.i.d. (different Vs)")
update_prior_covariance_rank1_ted_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_rank1_ted_notiid_rcpp <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_unconstrained_ted_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_ted_notiid(X,U,V,p)
update_prior_covariance_unconstrained_ted_notiid_rcpp <- function (X,U,V,p,...)
  update_prior_covariance_ted_notiid(X,U,V,p)

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al (2013).
# (This update can also be viewed as computing maximum-likelihood
# estimates of the parameters in the probabilistic PCA model.) Input p
# is a vector of weights associated with the rows of X. Input r
# specifies an optional constraint on U; when r < n, where U is an n x
# n matrix, at most r of the eigenvalues are positive in the updated
# matrix.
ted <- function (X, p, minval = 0, r = ncol(X)) {
    
  # Account for the weights by scaling the rows of X.
  T <- crossprod(sqrt(safenormalize(p))*X)
  
  # Find U maximizing the likelihood under the model x ~ N(0,U + I),
  # subject to U being a symmetric positive definite (s.p.d.) matrix,
  # and with the additional constraint that at most r of its
  # eigenvalues are positive.
  return(shrink_cov(T,minval,r))
}
