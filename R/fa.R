# Perform an M-step update for unconstrained prior covariance (U) using factor
# analysis for the special case when V = I for all samples; 
# that is, the model is x ~ N(0,U + I).
update_prior_covariance_unconstrained_fa_iid <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U, fa_unconstrained(X, U$mat, p))

# Perform an M-step update for unconstrained prior covariance (U) using factor
# analysis, allowing for different residual variances V among the samples.
update_prior_covariance_unconstrained_fa_notiid <- function (X, U, V, p, ...)
  stop("unconstrained_fa_notiid is not yet implemented")

# Perform an M-step update for a scaled prior covariance matrix (U) using 
# factor analysis for the special case when V = I for all samples; 
# that is, the model is x ~ N(0,sU + I).
update_prior_covariance_scaled_fa_iid <- function (X, U, p, ...){
  r <- sum(eigen(U$U0)$values > 1e-15)
  update_prior_covariance_struct_scaled(U, fa_scaled_iid(X, U$U0, p, U$s))
}

# Perform an M-step update for a scaled prior covariance (U) using factor
# analysis, allowing for different residual variances V among the samples.
update_prior_covariance_scaled_fa_notiid <- function (X, U, V, p, ...){
  r <- sum(eigen(U$U0)$values > 1e-15)
  update_prior_covariance_struct_scaled(U, fa_scaled_notiid(X, U$U0, V, p, U$s))
}

# Perform an M-step update for a rank1 prior covariance matrix U using
# factor analysis for the special case when V = I for all samples;
# that is, the model is x ~ N(0,U + I).
update_prior_covariance_rank1_fa_iid <- function (X, U, p, ...)
  update_prior_covariance_struct_rank1(U, fa_rank1_iid(X, U$vec, p))


# Perform an M-step update for a rank1 prior covariance matrix U using
# factor analysis allowing for different residual variances V among
# the samples.
update_prior_covariance_rank1_fa_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_struct_rank1(U, fa_rank1_notiid(X, U$vec, V, p))
  
# Perform an M-step update for unconstrained prior covariance matrix U
# using factor analyzer in the special case of V_j = I.  (See eq.(74)
# in main writeup.)
# @param X contains observed data of size n \times r.
# @param U the estimate of U in previous iteration
# @param p is a vector of weights
fa_unconstrained <- function(X, U, p){

  n <- nrow(X)
  r <- nrow(U)
  Q <- get_mat_Q(U)
  r <- ncol(Q)
  I <- diag(r)
  
  Sigma <- solve(crossprod(Q) + I)
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
                     
# Perform an M-step update for estimating the scalar for prior
# covariance matrix U0 in the special case of V_j = I. U0 can be
# rank-deficient.
# @param X contains observed data of size n \times r.
# @param U0 A known canonical covariance of size r \times r. 
# @param p is a vector of weights
# @param s is the scalar estimate in previous iteration
# @param r is the rank of U0
fa_scaled_iid<- function(X, U0, p, s){

  n = nrow(X)
  Q = get_mat_Q(U0)
  r = ncol(Q)
  I = diag(r)
  Sigma = solve(crossprod(Q) + I/s)
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

# Perform an M-step update for estimating the scalar for prior
# covariance matrix U0 in the general case where V_j can vary for
# different observations. U0 can be rank-deficient.
# @param X contains observed data of size n \times r.
# @param U0 A known canonical covariance of size r \times r. 
# @param V is a 3-d array, in which V[,,j] is the covariance matrix
#   for the jth observation
# @param p is a vector of weights
# @param s is the scalar estimate in previous iteration
# @param r is the rank of U0
fa_scaled_notiid <- function(X, U0, V, p, s){
  
  n = nrow(X)
  Q = get_mat_Q(U0)
  r = ncol(Q)
  I = diag(r)

  Sigma = c()
  V.inverse = c() 
  VinvQ = c()
  b = c()
  B = c()
  trB = rep(NA, n)  # trace of Bmat 
  
  for (i in 1:n){
    V.inverse[[i]] = solve(V[,,i])
    VinvQ[[i]] = solve(V[,,i]) %*% Q

    Sigma[[i]] = solve(crossprod(Q, VinvQ[[i]]) + I/s)    #t(Q) %*% VinvQ[[i]])
    b[[i]] = crossprod(X[i, ], VinvQ[[i]]) %*% Sigma[[i]]
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

# Perform an M-step update for rank1 prior covariance matrix U 
# using factor analyzer in the special case of V_j = I. 
# @param X contains observed data of size n \times r.
# @param u the estimate of vector that forms the rank1 U in previous iteration
# @param p is a vector of weights
fa_rank1_iid <- function (X, u, p) {
  n      <- nrow(X)
  m      <- ncol(X)
  sigma2 <- drop(1/(crossprod(u) + 1))
  mu <- sigma2*drop(tcrossprod(u, X)) # length-n vector
  theta <- mu*p
  eta <- p*(mu^2+sigma2)
  u = 1/sum(eta)* colSums(theta*X)
  return(u)
}

# Perform an M-step update for a prior rank-1 matrix, allowing for
# residual covariance matrices that differ among the data samples
# (rows of X).
# @param p is a vector of weights
# @param u is a vector
# @param V is a 3-d array, in which V[,,j] is the covariance matrix
# for the jth observation
fa_rank1_notiid<- function (X, u, V, p) {
  n      <- nrow(X)
  m      <- ncol(X)
  sigma2 <- rep(0,n)
  mu     <- rep(0,n)
  uw     <- matrix(0,m,n)
  Vinvw  <- array(0,c(m,m,n))
  for (i in 1:n){
    A          <- solve(V[,,i])
    utA        <- crossprod(u, A)
    sigma2     <- drop(1/(utA %*% u + 1))    # t(u) %*% A
    mu         <- drop(sigma2 * utA %*% X[i,])
    uw[,i]     <- p[i]*mu*A %*% X[i,]
    Vinvw[,,i] <- p[i]*(mu^2 + sigma2)*A
  }
  return(drop(solve(sliceSums(Vinvw),rowSums(uw))))
}
         



