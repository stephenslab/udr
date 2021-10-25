# Update the data structure for an unconstrained prior covariance
# matrix. Input "U" is the current data structure, and "mat" is the
# newly estimated matrix. This function is called in the
# update_prior_covariance_unconstrained_* functions below.
update_prior_covariance_unconstrained_struct <- function (U, mat) {
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$mat <- mat
  return(U)
}

# Update the data structure for a rank-1 prior covariance matrix.
# Input U is the current data structure, and "vec" is a vector
# containing the new estimates, such that the new rank-1 matrix is
# tcrossprod(vec). This function is called in the
# update_prior_covariance_rank1_* functions below.
update_prior_covariance_rank1_struct <- function (U, vec) {
  mat <- tcrossprod(vec)
  names(vec) <- names(U$vec)
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$vec <- vec
  U$mat <- mat
  return(U)
}

# Update the data structure for a scaled prior covariance matrix.
# Input U is the current data structure, and "s" is the new estimate
# of the scaling factor. This function is called in the
# update_prior_covariance_scaled_* functions below.
update_prior_covariance_scaled_struct <- function (U, s) {
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

# Perform an M-step update for a prior covariance matrix U using the
# "eigenvalue truncation" technique described in Won et al (2013).
# Note that input U is not used, and is included only for consistency
# with the other update_prior_covariance functions. Input p is a
# vector of weights associated with the rows of X. Input r specifies
# an optional constraint on U; when r < n, where U is an n x n matrix,
# at most r of the eigenvalues are positive in the updated matrix.
update_prior_covariance_unconstrained_ted <- function (X, U, V, p, minval,
                                                       r = ncol(X)) {
  if (!is.matrix(V))
    stop("unconstrained.update = \"ted\" does not work for case when data ",
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
  return(update_prior_covariance_unconstrained_struct(U,mat))
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ted.
update_prior_covariance_unconstrained_ted_rcpp <- function (X, U, V, p,
                                                             minval) {
  if (!is.matrix(V))
    stop("unconstrained.update = \"ted\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
  mat <- update_prior_covariance_ted_iid_rcpp(X,U$mat,V,p,minval)
  return(update_prior_covariance_unconstrained_struct(U,mat))
}

                     
#' Perform an M-step update for unconstrained prior covariance matrix U 
#' using factor analyzer in the special case of V_j = I. 
#' (See eq.(74) in main write-up)
#' @param X contains observed data of size n \times r.
#' @param U the estimate of U in previous iteration
#' @param p is a vector of weights
update_prior_covariance_unconstrained_fa <- function(X, U, p){
  n = nrow(X)
  r = nrow(U)
  Q = get_mat_Q(U, r)
  I = diag(r)
  
  Sigma = solve(t(Q) %*% Q + I)
  bmat = X %*% Q %*% Sigma   # n by r matrix
  
  A = t(X) %*% (p*bmat)
  B =  t(bmat) %*% (p*bmat) + sum(p)*Sigma
  Q = A %*% solve(B)
  return(Q%*% t(Q))
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
update_prior_covariance_rank1_fa <- function (X, U, V, p, minval) { 
  if (is.matrix(V)) {
    n <- nrow(X)
    m <- ncol(X)
    V <- array(V,c(m,m,n))
    vec <- update_prior_covariance_rank1_fa_iid(X,U$vec,p)
  } else 
    vec <- update_prior_covariance_rank1_fa_general(X,U$vec,V,p)
  return(update_prior_covariance_rank1_struct(U,vec))
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
# update_prior_covariance_unconstrained_ted for more information
# about the inputs.
update_prior_covariance_rank1_ted <- function (X, U, V, p, minval) {
  if (!is.matrix(V))
    stop("rank1.update = \"ted\" does not work for case when data ",
         "points are not i.i.d. (different Vs)")
  mat <- update_prior_covariance_unconstrained_ted(X,U,V,p,minval,r = 1)$mat
  vec <- getrank1(mat)
  return(update_prior_covariance_rank1_struct(U,vec))
}

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ted.
update_prior_covariance_rank1_ted_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_rank1_ted_rcpp is not yet implemented")
}




# Perform an M-step update for a prior rank-1 matrix, allowing for
# residual covariance matrices that differ among the data samples
# (rows of X).
# @param p is a vector of weights
# @param u is a vector
# @param V is a 3-d array, in which V[,,j] is the covariance matrix
# for the jth observation
update_prior_covariance_rank1_fa_general <- function (X, u, V, p) {
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
                   
                     
#' Perform an M-step update for rank1 prior covariance matrix U 
#' using factor analyzer in the special case of V_j = I. 
#' @param X contains observed data of size n \times r.
#' @param u the estimate of vector that forms the rank1 U in previous iteration
#' @param p is a vector of weights
update_prior_covariance_rank1_fa_iid <- function (X, u, p) {
  n      <- nrow(X)
  m      <- ncol(X)
  sigma2 <- drop(1/(t(u) %*% u + 1))
  mu <- sigma2*drop(t(u) %*% t(X)) # length-n vector
  theta <- mu*p
  eta <- p*(mu^2+sigma2)
  
  u = 1/sum(eta)* colSums(theta*X)
  return(u)
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
update_prior_covariance_scaled_fa <- function (X, U, V, p, minval) {
  if (is.matrix(V)){
    s <- update_prior_covariance_scaled_iid(X,U$U0,V,p,minval)
    }
  else{
    r <- sum(eigen(U$U0)$values > 1e-15)
    s <- update_prior_covariance_scaled_fa_general(X, U$U0, V, p, U$s, r)
    }
  return(update_prior_covariance_scaled_struct(U,s))
}

# This is a more efficient C++ implementation of
# update_prior_covariance_scaled_fa.
update_prior_covariance_scaled_fa_rcpp <- function (X, U, V, p, minval) {
  stop("update_prior_covariance_scaled_fa_rcpp is not yet implemented")
}

# Update the scaling factor for a prior canonical covariance (U) matrix.
# @param U0 A (fixed) covariance matrix.
# @param V A covariance matrix.
# @return An integer scalar
#
#' @importFrom stats uniroot
update_prior_covariance_scaled_iid <- function (X, U0, V, p, minval) {

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
                 c(0,1e6))$root)
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
                     
#' Perform an M-step update for estimating the scalar for prior covariance matrix U0
#' in the general case where V_j can vary for different observations. U0 can be rank-deficient.
#' @param X contains observed data of size n \times r.
#' @param U0 A known canonical covariance of size r \times r. 
#' @param V is a 3-d array, in which V[,,j] is the covariance matrix
# for the jth observation
#' @param p is a vector of weights
#' @param s is the scalar estimate in previous iteration
#' @param r is the rank of U0
update_prior_covariance_scaled_fa_general <- function(X, U0, V, p, s, r){
  
  n = nrow(X)
  Q = get_mat_Q(U0, r)
  I = diag(r)
  
  Sigma = c()
  V.inverse = c()
  b = c()
  B = c()
  trB = rep(NA, n)  # trace of Bmat 
  
  for (i in 1:n){
    V.inverse[[i]] = solve(V[,,i])
    Sigma[[i]] = solve(t(Q)%*% V.inverse[[i]] %*% Q + I/s)
    b[[i]] = t(t(X[i, ])%*% V.inverse[[i]] %*% Q %*% Sigma[[i]])
    B[[i]] = b[[i]] %*% t(b[[i]])+Sigma[[i]]
    trB[i] = sum(diag(B[[i]]))
  }
  
  if (r ==1){
    s = sum(p*B)/sum(p)
  }else{
    s = sum(trB*p)/(sum(p)*r)
  }
  
  return(s)
}
                     
                     
#' Perform an M-step update for estimating the scalar for prior covariance matrix U0
#' in the special case of V_j = I. U0 can be rank-deficient.
#' @param X contains observed data of size n \times r.
#' @param U0 A known canonical covariance of size r \times r. 
#' @param p is a vector of weights
#' @param s is the scalar estimate in previous iteration
#' @param r is the rank of U0
update_prior_covariance_scaled_fa_iid <- function(X, U0, p, s, r){
  
  n = nrow(X)
  Q = get_mat_Q(U0, r)
  I = diag(r)
  
  Sigma = solve(t(Q) %*% Q + I/s)
  bmat = X %*% Q %*% Sigma   # n by r matrix to store b_j
  B = c()
  trB = rep(NA, n)  # trace of Bmat 
  
  for (i in 1:n){
    B[[i]] = bmat[i, ] %*% t(bmat[i, ])+ Sigma
    trB[i] = sum(diag(B[[i]]))
  }
  
  if (r ==1){
    s = sum(p*B)/sum(p)
  }else{
    s = sum(trB*p)/(sum(p)*r)
  }
  return(s)
}
