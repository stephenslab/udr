#' Perform an M-step update for unconstrained prior covariance matrix U 
#' using factor analyzer in the special case of V_j = I. 
#' (See eq.(74) in main write-up)
#' @param X contains observed data of size n \times r.
#' @param U the estimate of U in previous iteration
#' @param p is a vector of weights
update_prior_covariance_unconstrained_fa <- function(X, U, p){
  n <- nrow(X)
  r <- nrow(U)
  Q <- get_mat_Q(U, r)
  I <- diag(r)
  
  Sigma <- solve(t(Q) %*% Q + I)
  bmat <- X %*% Q %*% Sigma   # n by r matrix
  A <- crossprod(X, p*bmat) # t(X) %*% (p*bmat)
  B <-  crossprod(bmat, p*bmat) + sum(p)*Sigma
  Q <- A %*% solve(B)
  mat <- tcrossprod(Q)
  return(mat)
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
  sigma2 <- drop(1/(crossprod(u) + 1))
  mu <- sigma2*drop(tcrossprod(u, X)) # length-n vector
  theta <- mu*p
  eta <- p*(mu^2+sigma2)
  u = 1/sum(eta)* colSums(theta*X)
  return(u)
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
    b[[i]] = t(X[i, ])%*% V.inverse[[i]] %*% Q %*% Sigma[[i]]
    B[[i]] = crossprod(b[[i]])+ Sigma[[i]]
    trB[i] = sum(diag(B[[i]]))
  }
  
  if (r ==1){
    s = sum(p*unlist(B))/sum(p)
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
    B[[i]] = tcrossprod(bmat[i, ]) + Sigma
    trB[i] = sum(diag(B[[i]]))
  }
  
  if (r ==1){
    s = sum(p*unlist(B))/sum(p)
  }else{
    s = sum(trB*p)/(sum(p)*r)
  }
  return(s)
}
