# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), for the special case when V =
# I for all samples; that is, the model is x ~ N(0,U + I).
update_prior_covariance_unconstrained_ed_iid <- function (X, U, p, n0, S0, ...){
  res <- ed_reg_iid(X, U$mat, p, S0, n0, U$s)
  return(update_prior_covariance_struct_unconstrained(U, res$U, res$sigma2))
}
  

# This is a more efficient C++ implementation of 
# update_prior_covariance_unconstrained_ed_iid.
update_prior_covariance_unconstrained_ed_iid_rcpp <- function (X, U, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed_iid_rcpp(X,U$mat,p))

# Perform an M-step update for a prior covariance (U) using the update
# formula derived in Bovy et al (2011), allow for different residual
# variances V among the samples.
update_prior_covariance_unconstrained_ed_notiid <- function (X, U, V, p, ...)
  update_prior_covariance_struct_unconstrained(U,ed(X,U$mat,V,p))

# This is a more efficient C++ implementation of
# update_prior_covariance_rank1_ed_notiid. (The C++ version has not yet been
# implemented, so for now we simply call the R implementation.)
update_prior_covariance_unconstrained_ed_notiid_rcpp <- function (X, U, V, p,
                                                                  ...) {
  message("update_prior_covariance_unconstrained_ed_notiid_rcpp is not yet ",
          "implemented; using R version instead")
  return(update_prior_covariance_unconstrained_ed_notiid(X,U,V,p,...))
}

# These functions are defined only to provide more informative error
# messages.
update_prior_covariance_ed_invalid <- function (X, U, V, p, ...)
  stop("control$scale.update = \"ed\" and control$rank1.update = \"ed\" ",
       "are not valid choices")
update_prior_covariance_scaled_ed_iid         <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_iid_rcpp    <- function (X, U, p,...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_iid          <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_iid_rcpp     <- function (X, U, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_notiid      <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_scaled_ed_notiid_rcpp <- function (X, U, V, p,...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_notiid       <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()
update_prior_covariance_rank1_ed_notiid_rcpp  <- function (X, U, V, p, ...)
  update_prior_covariance_ed_invalid()

# Update the prior covariance matrix (U) in the model x ~ N(0,U + I)
# using the update formula derived in Bovy et al (2011). Input p is a
# vector of weights associated with the rows of X.
ed_iid <- function (X, U, p) {
  m  <- ncol(X)
  I  <- diag(m)
  B  <- solve(U + I,U)
  X1 <- crossprod((sqrt(safenormalize(p))*X) %*% B)
  return(U - U %*% B + X1)
}


#' Function for scaled regularized ED with an inverse Wishart
#' prior on U.tilde. U.tilde ~ W_R^{-1}(S0, nu), where nu = n0-R-1. 
#' U = sigma^2*U.tilde.
#' The derivation is available in write-up: Notes on estimation of 
#' large covariance matrices.
#' @param X: data matrix of size $n$ by $R$.
#' @param U: initialization of U or estimate from previous iteration
#' @param p: the weight vector for a component
#' @param S0: a covariance matrix in inverse-Wishart prior
#' @param n0: parameter in inverse-Wishart prior, n0 = nu + R + 1
#' @param sigma2: initialization of the scalar value or estimate from previous iteration.
ed_reg_iid <- function(X, U, p, S0, n0, sigma2){
  m  <- ncol(X)
  I  <- diag(m)
  A  <- solve(U + I,U)
  B <- U %*% (I-A)
  bmat <- X %*% A
  U <- (crossprod(sqrt(p)*bmat)+sum(p)*B+sigma2*n0*S0)/(sum(p)+ n0)
  if (n0 != 0)
    sigma2 <- m/sum(diag(S0 %*% solve(U)))
  return(list(U = U, sigma2 = sigma2))
}

# Perform an M-step update for a prior covariance matrix (U) using the
# update formula derived in Bovy et al (2011), allowing for residual
# covariance matrices V that differ among the data samples (rows of X).
ed <- function (X, U, V, p) {
  n  <- nrow(X)
  m  <- ncol(X)
  B  <- array(0,c(m,m,n))
  bb <- array(0,c(m,m,n))
  for (i in 1:n) {
    Tinv    <- solve(U + V[,,i])
    bb[,,i] <- tcrossprod(((sqrt(p[i])*U) %*% Tinv) %*% X[i,])
    B[,,i]  <- p[i]*(U - (U %*% Tinv) %*% U)
  }
  return((sliceSums(bb) + sliceSums(B))/sum(p))
}
