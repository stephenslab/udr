# Perform an M-step update for the mixture weights in the mixture
# prior (or no update if update = "none")
update_mixture_weights_em <- function (P, update) {
  if (update == "em")
    wnew <- colSums(P)/nrow(P)
  else if (update == "none")
    wnew <- w
  else
    stop("control$weights.update == \"",update,"\" is not implemented")
  return(wnew)
}

# Perform an M-step update for the residual covariance matrix (or no
# update if update = "none"). Input argument P should be the n x k
# matrix of posterior mixture assignment probabilities returned by
# compute_posterior_probs. Input argument V should be an m x m
# matrix. Input argument U may either be a list of length k in which
# U[[i]]$mat is an m x m matrix, or an m x m x k array.
update_resid_covariance <- function (X, U, V, P, update,
                                     version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (is.list(U))
    U <- ulist2array(U)
  if (update == "em") {
    if (version == "R")
      Vnew <- update_resid_covariance_helper(X,U,V,P)
    else if (version == "Rcpp")
      Vnew <- update_resid_covariance_rcpp(X,U,V,P)
  } else if (update == "none")
    Vnew <- V
  else
    stop("control$resid.update == \"",update,"\" is not implemented")
  return(Vnew)
}

# Implements update_resid_covariance with version = "R".
update_resid_covariance_helper <- function (X, U, V, P) {
  n <- nrow(X)
  m <- ncol(X)
  Vnew <- matrix(0,m,m)
  for (i in 1:n) {
    out <- compute_posterior_mvtnorm_mix(X[i,],P[i,],U,V)
    Vnew <- Vnew + out$S1 + tcrossprod(X[i,] - out$mu1)
  }
  return(Vnew/n)
}

# Determine the functions for updating the prior covariance matrices
# based on the (1) the covariance matrix types and (2) the control
# settings. The return value is a character vector with one entry for
# each prior covariance matrix.
assign_prior_covariance_updates <- function (covtypes, control) {
  k <- length(covtypes)
  covupdates <- rep(as.character(NA),k)
  for (i in 1:k)
    covupdates[i] <- paste0("update_prior_covariance_",covtypes[i],"_",
                            control[[paste(covtypes[i],"update",sep = ".")]],
                            ifelse(control$version == "Rcpp","_rcpp",""))
  names(covupdates) <- names(covtypes)
  return(covupdates)
}

# Perform M-step updates for all the prior covariance matrices U.
# Input argument V may either be an m x m matrix, a list of m x m
# matrices of length n, or a m x m x n array.
update_prior_covariances <- function (X, U, V, P, covupdates, minval) {
  k <- length(U)
  if (is.list(V))
    V <- list2array(V)
  for (i in 1:k)
    U[[i]] <- do.call(covupdates[i],list(X = X,U = U[[i]],V = V,p = P[,i],
                                         minval = minval))
  return(U)
}

# This function simply returns the unconstrained prior covariance
# matrix without updating it.
update_prior_covariance_unconstrained_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the unconstrained prior covariance
# matrix without updating it.
update_prior_covariance_unconstrained_none_rcpp <- function (X, U, V, p,
                                                             minval) {
  return(U)
}

# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011). Input p is a vector of weights
# associated with the rows of X.
update_prior_covariance_unconstrained_ed <- function (X, U, V, p, minval) {
  if (is.matrix(V))
    U$mat <- update_prior_covariance_ed_iid(X,U$mat,V,p)
  else
    U$mat <- update_prior_covariance_ed_general(X,U$mat,V,p)
  return(U)
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed.
update_prior_covariance_unconstrained_ed_rcpp <- function (X, U, V, p,
                                                           minval) {
  if (is.matrix(V)) 
    U$mat <- update_prior_covariance_ed_iid_rcpp(X,U$mat,V,p)
  else
    stop("update_prior_covariance_unconstrained_ed_rcpp is not yet ",
         "implemented for case when data points are not i.i.d. (different Vs)")
  return(U)
}

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al (2013).
# Note that input U is not used, and is included only for consistency
# with the other update_prior_covariance functions. Input p is a
# vector of weights associated with the rows of X. Input r specifies
# an optional constraint on U; when r < n, where U is an n x n matrix,
# at most r of the eigenvalues are positive in the updated matrix.
update_prior_covariance_unconstrained_teem <- function (X, U, V, p, minval,
                                                        r = ncol(X)) {
  if (!is.matrix(V))
    stop("unconstrained.update = \"teem\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
    
  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  p <- safenormalize(p)
  R <- chol(V)
  T <- crossprod((sqrt(p)*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite, with at most r of its eigenvalues being
  # positive.
  U$mat <- shrink_cov(T,minval,r)
  
  # Recover the solution for the original (untransformed) data.
  U$mat <- t(R) %*% U$mat %*% R
  return(U)
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_teem.
update_prior_covariance_unconstrained_teem_rcpp <- function (X, U, V, p,
                                                             minval) {
  if (!is.matrix(V))
    stop("unconstrained.update = \"teem\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
  U$mat <- update_prior_covariance_teem_iid_rcpp(X,U$mat,V,p,minval)
  return(U)
}

# This function simply returns the rank-1 prior covariance matrix
# without updating it.
update_prior_covariance_rank1_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the rank-1 prior covariance matrix
# without updating it.
update_prior_covariance_rank1_none_rcpp <- function (X, U, V, p, minval) {
  return(U)
}

# Perform an M-step update for a rank-1 prior covariance (U) using the
# update formula derived originally by David Gerard. Input p is a
# vector of weights associated with the rows of X.
update_prior_covariance_rank1_ed <- function (X, U, V, p, minval) { 
  if (is.matrix(V)) {
    n <- nrow(X)
    m <- ncol(X)
    V <- array(V,c(m,m,n))
    U$vec <- update_prior_covariance_rank1_ed_general(X,U$vec,V,p)
  } else 
    U$vec <- update_prior_covariance_rank1_ed_general(X,U$vec,V,p)
  U$mat <- tcrossprod(U$vec)
  return(U)
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ed.
update_prior_covariance_rank1_ed_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_rank1_ed_rcpp is not yet implemented") 
}

# Perform an M-step update for a rank-1 prior covariance matrix U
# using the "eigenvalue truncation" technique. Note that input U is
# not used, and is included only for consistency with the other
# update_prior_covariance functions. See
# update_prior_covariance_unconstrained_teem for more information
# about the inputs.
update_prior_covariance_rank1_teem <- function (X, U, V, p, minval) {
  if (!is.matrix(V))
    stop("rank1.update = \"teem\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
  U <- update_prior_covariance_unconstrained_teem(X,U,V,p,minval,r = 1)
  U$vec <- getrank1(U$mat)
  return(U)
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_teem.
update_prior_covariance_rank1_teem_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_rank1_teem_rcpp is not yet implemented")
}

# This function simply returns the scaled prior covariance matrix
# without updating it.
update_prior_covariance_scaled_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the scaled prior covariance matrix
# without updating it.
update_prior_covariance_scaled_none_rcpp <- function (X, U, V, p, minval) {
  return(U)
}

# Perform an M-step update for a scaled prior ocvariance matrix (U).
# Input p is a vector of weights associated with the rows of X.
update_prior_covariance_scaled_em <- function (X, U, V, p, minval) {
  if (is.matrix(V))
    U$s <- update_prior_covariance_scaled_em_iid(X,U$U0,V,p,minval)
  else 
    stop("update_prior_covariance_scaled_em is not yet implemented for ",
         "case when V is an array")
  U$mat <- with(U,s*U0)
  return(U)
}

# This is a more efficient C++ implementation of
# update_prior_covariance_scaled_em.
update_prior_covariance_scaled_em_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_scaled_em_rcpp is not yet implemented")
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
# when the residual covariances V are the same for all data
# points. Input p is a vector of weights associated with the rows of
# X.
update_prior_covariance_ed_iid <- function (X, U, V, p) {
  p <- safenormalize(p)
  T <- U + V
  B <- solve(T,U)
  X1 <- crossprod((sqrt(p)*X) %*% B)
  return(U - U %*% B + X1)
}

# Perform an M-step update for a prior rank-1 matrix, allowing for
# residual covariance matrices that differ among the data samples
# (rows of X).
# @param p is a vector of weights
# @param u is a vector
# @param V is a 3-d array, in which V[,,j] is the covariance matrix
# for the jth observation
update_prior_covariance_rank1_ed_general <- function (X, u, V, p) {
  n      <- nrow(X)
  m      <- ncol(X)
  sigma2 <- rep(0,n)
  mu     <- rep(0,n)
  uw     <- matrix(0,m,n)
  Vinvw  <- array(0,c(m,m,n))
  for (i in 1:n){
    A          <- solve(V[,,i])
    sigma2     <- drop(1/(t(u) %*% A %*% u + 1))
    mu         <- drop(sigma2*t(u) %*% A %*% X[i,])
    Vinvw[,,i] <- p[i]*(mu^2 + sigma2)*A
    uw[,i]     <- p[i]*mu*A %*% X[i,]
  }
  return(drop(solve(sliceSums(Vinvw),rowSums(uw))))
}

# Update the scaling factor for a prior canonical covariance (U) matrix.
# @param U0 A (fixed) covariance matrix.
# @param V A covariance matrix.
# @return An integer scalar
#
#' @importFrom stats uniroot
update_prior_covariance_scaled_em_iid <- function (X, U0, V, p, minval) {

    # Transform data using the trick
    # Uhat = R^{-T}UR^{-1}
    # Xhat = XR^{-1}
    R    <- chol(V)
    Uhat <- t(solve(R))%*% U0 %*% solve(R)
    Xhat <- X %*% solve(R)
    
    # Eigenvalue decomposition based on transformed U.
    evd <- eigen(Uhat)
    lambdas <- ifelse(evd$values < minval,0,evd$values)

    Y <- t(evd$vectors) %*% t(Xhat)  # Y: p by n
    return(uniroot(function(s) grad_loglik_scale_factor(s,p,Y,lambdas),c(0,1),
                   extendInt = "yes")$root)
}

# Function for 1-d search of s value based on eq. (20) in the write-up
# @param p Weights for one component
# @param s The scaling factor for one component we aim to search for
# @param Y The transformed data
# @param lambdas Eigenvalues of U.
# @return A function of the scalar s.
grad_loglik_scale_factor <- function (s, p, Y, lambdas) {
  unweighted_sum <-
    apply(Y,2,function(y) sum(lambdas*y^2/((s*lambdas + 1)^2)) -
                          sum(lambdas/(s*lambdas + 1)))
  return(sum(p*unweighted_sum))
}
