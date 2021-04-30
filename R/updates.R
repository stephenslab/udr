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
# update if update = "none").
update_resid_covariance <- function (X, U, V, P, update,
                                     version = c("Rcpp","R")) {
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

# Perform an M-step update for the prior covariance matrices.
update_prior_covariances <- function (X, U, V, P, covtypes, control) {
  if (control$version == "Rcpp") {
    control$version <- "R"
    message("update_prior_covariances with version = \"Rcpp\" is not yet ",
            "implemented; switching to version = \"R\"")
  }
  if (control$version == "R")
    Unew <- update_prior_covariances_helper(X,U,V,P,covtypes,control)
  return(Unew)
}

# Implements update_prior_covariances with version = "R".
update_prior_covariances_helper <- function (X, U, V, P, covtypes,
                                             control) {
  Unew <- U
  k <- ncol(P)
  for (i in 1:k) {
    if (covtypes[i] == "scaled") {

      # Update the scaling factor.
      if (control$scaled.update == "em") {
        if (!is.matrix(V))
          stop("control$scaled.update == \"em\" can only be used when ",
               "the residual covariance (V) is the same for all data points")
          Unew[,,i] <- U[,,i] *
            update_prior_scalar(X,U[,,i],V,P[,i],control$minval)
      } else if(control$scaled.update != "none")
        stop("control$scaled.update == \"",control$scaled.update,
             "\" is not implemented")
    } else if (covtypes[i] == "rank1") {

      # Update a rank-1 covariance matrix.
    if (control$rank1.update == "teem") {
      if (!is.matrix(V))
        stop("control$rank1.update == \"teem\" is currently not ",
             "implemented for case when each data point has a different ",
             "residual covariance, V")
      Unew[,,i] <- update_prior_covariance_teem(X,V,P[,i],control$minval, r = 1)
    } else if (control$rank1.update == "ed"){
        if (!is.matrix(V)){
            Unew[,,i] <- update_prior_rank1_general(X,U[,,i],V,P[,i])
        }
    }else if (control$rank1.update != "none"){
        stop("control$rank1.update == \"",control$rank1.update,
        "\" is not implemented")
    }
           
    }else if (covtypes[i] == "unconstrained") {

      # Update a full (unconstrained) matrix.
      if (control$unconstrained.update == "ed") {
        if (!is.matrix(V))
          Unew[,,i] <- update_prior_covariance_ed_general(X,U[,,i],V,P[,i])
        else
          Unew[,,i] <- update_prior_covariance_ed(X,U[,,i],V,P[,i])
      } else if (control$unconstrained.update == "teem") {
        if (!is.matrix(V))
          stop("control$unconstrained.update == \"teem\" can only be used ",
               "when the residual covariance (V) is the same for all data ",
               "points")
        Unew[,,i] <- update_prior_covariance_teem(X,V,P[,i],control$minval, r = nrow(V))
      } else if (control$unconstrained.update != "none")
        stop("control$unconstrained.update == \"",control$unconstrained.update,
             "\" is not implemented")
    } else
      stop("Invalid prior covariance type")
  }
  return(Unew)
}

# Perform an M-step update for one of the prior covariance matrices
# using the update formula derived in Bovy et al (2011). Here, p is a
# vector, with one entry per row of X, giving the posterior assignment
# probabilities for the mixture component being updated.
update_prior_covariance_ed <- function (X, U, V, p) {
  p <- safenormalize(p)
  T <- U + V
  B <- solve(T,U)
  X1 <- crossprod((sqrt(p) * X) %*% B)
  return(U - U %*% B + X1)
}

# Perform an M-step update for one of the prior covariance matrices
# using the eigenvalue-truncation technique described in Won et al
# (2013). Input p is a vector, with one entry per row of X, giving
# the posterior assignment probabilities for the mixture components
# being updated.
update_prior_covariance_teem <-
  function (X, V, p, minval, r) {

  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  p <- safenormalize(p)
  R <- chol(V)
  T <- crossprod((sqrt(p)*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite, or U being a rank-1 matrix. The unconstrained
  # update is based on the fact that the covariance matrix that minimizes the
  # likelihood subject to the constraint that U is positive definite is obtained by
  # truncating the eigenvalues of T = U + I less than 1 to be 1; see
  # Won et al (2013), p. 434, the sentence just after equation (16).
  # For rank1 update, only the first rth largest eigenvalues of U are kept.
  # The others are set to 0.

  U <- shrink_cov(T, minval, r)
    
  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}

# "Shrink" matrix T = U + I for (1) unconstrained update and (2) rank1 update.
# For unconstrained update, find the "best" matrix T
# satisfying the constraint that U is positive definite (if minval >
# 0) or positive semi-definite (if minval <= 0). This is achieved by
# setting any eigenvalues of T less than 1 to 1 + minval, or,
# equivalently, setting any eigenvalues of U less than 0 to be
# minval. The output is a positive definite matrix, U, or a positive
# semi-definite matrix if minval <= 0.
# For rank1 update, we keep only first rth eigenvalues where r is the rank of matrix U.
shrink_cov <- function (T, minval = 0, r) {
    
  minval <- max(0,minval)
  out <- eigen(T)
  d <- out$values
  
  if (nrow(T)!= r){
      d[(r+1):length(d)] <- 1
  }
  d <- pmax(d - 1,minval)
  return(tcrossprod(out$vectors %*% diag(sqrt(d))))
}

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
    return(uniroot(function(s) grad_loglik_scale_factor(s,p,Y,lambdas),c(0,1e6),
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
