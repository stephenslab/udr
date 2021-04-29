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
update_prior_covariances <- function (X, U, V, P, covupdates, minval) {
  k <- length(U)
  for (i in 1:k)
    U[[i]] <- do.call(covupdates[i],list(X = X,U = U[[i]],V = V,p = P[,i],
                                         minval = minval))
  return(U)
}

# This function simply returns the scaled prior covariance matrix
# without updating it.
update_prior_covariance_scaled_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the rank-1 prior covariance matrix
# without updating it.
update_prior_covariance_rank1_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the unconstrained prior covariance
# matrix without updating it.
update_prior_covariance_unconstrained_none <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the scaled prior covariance matrix
# without updating it.
update_prior_covariance_scaled_none_rcpp <- function (X, U, V, p, minval) {
  return(U)
}

# This function simply returns the rank-1 prior covariance matrix
# without updating it.
update_prior_covariance_rank1_none_rcpp <- function (X, U, V, p, minval) {
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
    stop("update_prior_covariance_unconstrained_ed is not yet implemented ",
         "for case when data points are not i.i.d. (different Vs)")
  return(U)
}

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed.
update_prior_covariance_unconstrained_ed_rcpp <- function (X, U, V, p,
                                                           minval) {
  if (is.matrix(V)) 
    U$mat <- update_prior_covariance_ed_iid_rcpp(X,U$mat,V,p)
  else
    stop("update_prior_covariance_unconstrained_ed is not yet implemented ",
         "for case when data points are not i.i.d. (different Vs)")
  return(U)
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

# update_prior_covariances_helper = function (X, U, V, P, covtypes, control) {
#     if (covtypes[i] == "scaled") {

#       # Update the scaling factor.
#       if (control$scaled.update == "em") {
#           Unew[,,i] <- U[,,i] *
#             update_prior_scalar(X,U[,,i],V,P[,i],control$minval)
#     } else if (covtypes[i] == "rank1") 

#       # Update a rank-1 covariance matrix.
#     if (control$rank1.update == "teem") {
#       Unew[,,i] <- update_prior_covariance_teem(X,V,P[,i],control$minval,r = 1)
#     } else if (control$rank1.update == "ed"){
#         if (!is.matrix(V)){
#             Unew[,,i] <- update_prior_rank1_general(X,U[,,i],V,P[,i])
#         }
        
#     } else if (covtypes[i] == "unconstrained") {

#       # Update a full (unconstrained) matrix.
#       if (control$unconstrained.update == "ed") {
#         if (!is.matrix(V))
#           Unew[,,i] <- update_prior_covariance_ed_general(X,U[,,i],V,P[,i])

# Update the scaling factor for prior canonical covariance matrix.
# @param U0 Canonical covariance matrix
# @return An integer scalar
update_prior_scalar <- function (X, U0, V, p, minval){

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

# Function for 1-d search of s value based on eq.(20) in the write-up
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

# Perform an M-step update for one of the prior covariance matrices
# using the update formula derived in Bovy et al (2011) with varied V_j.
# @param p is a vector of the weight matrix for one component.
# @param U is a matrix
# @param V is a 3-d array object, containing V_j for each observation
update_prior_covariance_ed_general = function(X, U, V, p){

    
  B.weighted = c()
  b.weighted = c()
  n <- nrow(X)
  
  for (i in 1:n){
    b.weighted[[i]] = sqrt(p[i])* U %*% solve(U+V[,,i]) %*% X[i, ]
    B.weighted[[i]] = p[i]*(U - U %*% solve(U+V[,,i])%*% U)
  }
  
  bb.weighted = lapply(b.weighted, function(x) as.matrix(x %*% t(x)))
  bb.weighted = simplify2array(bb.weighted)
  B.weighted = simplify2array(B.weighted)
  
  Unew = (apply(bb.weighted, c(1,2), sum) + apply(B.weighted, c(1,2), sum))/sum(p)
  return(Unew)
}

# Perform an M-step update for one of the prior rank1 matrix when V varies
# The algorithm is based on Bovy et al (2011) and derived by David Gerard.
# @param p is a vector of the weight matrix for one component.
# @param U is a matrix
# @param V is a 3-d array object, containing V_j for each observation
update_prior_rank1_general = function(X, U, V, p){
  
  n = nrow(X)
  vec = get_vec(U)
  
  sigma2 = rep(NA, n)
  mu = rep(NA, n)
  V.inverse = c()

  for (i in 1:n){
    V.inverse[[i]] = solve(V[,,i])
    sigma2[i] = 1/ (t(vec) %*%  V.inverse[[i]] %*% vec + 1)
    mu[i] = sigma2[i] * t(vec)%*% V.inverse[[i]] %*% X[i,]
  }
  
  theta = p*mu
  eta = p*(mu^2 + sigma2)
  
  vec.weighted = lapply(1:n, function(i) theta[i]*V.inverse[[i]]%*% X[i, ])
  V.inverse.weighted  = lapply(1:n, function(i) eta[i]*V.inverse[[i]])
  
  sum1 = apply(simplify2array(vec.weighted), c(1,2), sum)
  sum2 = apply(simplify2array(V.inverse.weighted), c(1,2), sum)
  
  vnew = solve(sum2) %*% sum1
  Unew = tcrossprod(vnew)
  
  return(Unew)
}

# function to get the vector from rank1 matrix.
get_vec = function(U){
  evd = eigen(U)
  vec =  sqrt(evd$values[1])*evd$vectors[, 1]
  return(vec)
}
