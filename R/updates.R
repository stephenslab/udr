# Perform an M-step update for the mixture weights in the mixture
# prior (or no update if update = "none")
#
#' @rdname ud_fit_advanced
#' 
#' @param update Describe input argument "update" here.
#' 
#' @export
#' 
update_mixture_weights_em <- function (fit, update = c("em","none")) {
  update <- match.arg(update)
    
  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\",",
         "such as the output of ud_init")
    
  # Update the mixture weights.
  if (update == "em")
    wnew <- colSums(fit$P)/nrow(fit$P)
  else if (update == "none")
    wnew <- fit$w
  else 
    stop("control$weights.update == \"",update,"\" is not implemented")

  # Add the names back.
  names(wnew) <- names(fit$w)
  
  # Output the updated fit.
  fit$w <- wnew
  return(fit)
}

# Perform an M-step update for the residual covariance matrix (or no
# update if update = "none"). Input argument P should be the n x k
# matrix of posterior mixture assignment probabilities returned by
# compute_posterior_probs. Input argument V should be an m x m
# matrix. Input argument U may either be a list of length k in which
# U[[i]]$mat is an m x m matrix, or an m x m x k array.
#
#' @rdname ud_fit_advanced
#' 
#' @export
#' 
update_resid_covariance <- function (fit, update = c("em","none"),
                                     version = c("Rcpp","R")) {
  update  <- match.arg(update)
  version <- match.arg(version)

  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\",",
         "such as the output of ud_init")
  if (!is.matrix(fit$V))
    stop("The residual covariance V can only be updated when it is the same ",
         "for all data points")

  # Process U.
  U <- ulist2array(fit$U)

  # Update the residual covariance matrix, V.
  if (update == "em") {
    if (version == "R")
      Vnew <- update_resid_covariance_helper(fit$X,U,fit$V,fit$P)
    else if (version == "Rcpp")
      Vnew <- update_resid_covariance_rcpp(fit$X,U,fit$V,fit$P)
  } else if (update == "none")
    Vnew <- fit$V
  else
    stop("control$resid.update == \"",update,"\" is not implemented")

  # Add back the row and column names.
  rownames(Vnew) <- rownames(fit$V)
  colnames(Vnew) <- colnames(fit$V)
  
  # Output the updated fit.
  fit$V <- Vnew
  return(fit)
}

# Perform M-step updates for all the prior covariance matrices U.
# Input argument V may either be an m x m matrix, a list of m x m
# matrices of length n, or a m x m x n array.
#
#' @rdname ud_fit_advanced
#'
#' @param minval Describe input argument "minval" here.
#' 
#' @export
#' 
update_prior_covariances <-
  function (fit, covupdates = rep("none",length(fit$U)), minval = 1e-14) {

  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\",",
         "such as the output of ud_init")

  # Process V.
  V <- fit$V
  if (is.list(V))
    V <- list2array(V)

  # Update the prior covariance matrices, U.
  k <- length(fit$U)
  for (i in 1:k) {
    fit$U[[i]] <- do.call(covupdates[i],list(X = fit$X,U = fit$U[[i]],V = V,
                                             p = fit$P[,i],minval = minval))
  }
  
  # Output the updated fit.
  return(fit)
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

# TO DO: Explain here what this function does, and how to use it.
update_prior_covariance_unconstrained <- function (U, mat) {
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$mat <- mat
  return(U)
}

# TO DO: Explain here what this function does, and how to use it.
update_prior_covariance_rank1 <- function (U, vec) {
  mat <- tcrossprod(vec)
  names(vec) <- names(U$vec)
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$vec <- vec
  U$mat <- mat
  return(U)
}

# TO DO: Explain here what this function does, and how to use it.
update_prior_covariance_scaled <- function (U, s) {
  U$s   <- s
  U$mat <- s * U$U0
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
    mat <- update_prior_covariance_ed_iid(X,U$mat,V,p)
  else
    mat <- update_prior_covariance_ed_general(X,U$mat,V,p)
  return(update_prior_covariance_unconstrained(U,mat))
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
  return(update_prior_covariance_unconstrained(U,mat))
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
  mat <- shrink_cov(T,minval,r)
  
  # Recover the solution for the original (untransformed) data.
  mat <- t(R) %*% mat %*% R
  return(update_prior_covariance_unconstrained(U,mat))
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_teem.
update_prior_covariance_unconstrained_teem_rcpp <- function (X, U, V, p,
                                                             minval) {
  if (!is.matrix(V))
    stop("unconstrained.update = \"teem\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
  mat <- update_prior_covariance_teem_iid_rcpp(X,U$mat,V,p,minval)
  return(update_prior_covariance_unconstrained(U,mat))
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
    vec <- update_prior_covariance_rank1_ed_general(X,U$vec,V,p)
  } else 
    vec <- update_prior_covariance_rank1_ed_general(X,U$vec,V,p)
  return(update_prior_covariance_rank1(U,vec))
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
  mat <- update_prior_covariance_unconstrained_teem(X,U,V,p,minval,r = 1)$mat
  vec <- getrank1(mat)
  return(update_prior_covariance_rank1(U,vec))
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

# Perform an M-step update for a scaled prior covariance matrix (U).
# Input p is a vector of weights associated with the rows of X.
update_prior_covariance_scaled_em <- function (X, U, V, p, minval) {
  if (is.matrix(V))
    s <- update_prior_covariance_scaled_em_iid(X,U$U0,V,p,minval)
  else 
    stop("update_prior_covariance_scaled_em is not yet implemented for ",
         "case when V is an array")
  return(update_prior_covariance_scaled(U,s))
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
# when the residual covariances V are the same for all data points.
# Input p is a vector of weights associated with the rows of X.
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
  # Uhat = R^{-T}*U*R^{-1}
  # Xhat = X*R^{-1}
  R    <- chol(V)
  Uhat <- solve(t(R),U0) %*% solve(R)
  Xhat <- X %*% solve(R)
    
  # Eigenvalue decomposition based on transformed U.
  evd <- eigen(Uhat)
  lambdas <- ifelse(evd$values < minval,minval,evd$values)

  Y <- t(Xhat %*% evd$vectors)
  return(uniroot(function (s) grad_loglik_scale_factor(s,p,Y,lambdas),
                 c(0,1),extendInt = "yes")$root)
}

# Function for 1-d search of s value based on eq. (20) in the write-up.
# @param p Vector of weights.
# @param s The scaling factor for one component we aim to search for
# @param Y The transformed data
# @param lambdas Eigenvalues of U.
# @return A function of the scalar s.
grad_loglik_scale_factor <- function (s, p, Y, lambdas)
  sum(p*apply(Y,2,function(y) sum(lambdas*y^2/((s*lambdas + 1)^2)) -
                              sum(lambdas/(s*lambdas + 1))))

