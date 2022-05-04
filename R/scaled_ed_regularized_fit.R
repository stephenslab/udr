
#' Function for scaled regularized ED with an inverse Wishart
#' prior on \tilde U. \tilde U \sim W_R^{-1}(S0, \nu), where \nu = n0-R-1. 
#' U = \sigma^2*\tilde U.
#' The derivation is available in write-up: Notes on estimation of 
#' large covariance matrices.
#' @param X: data matrix of size $n$ by $R$.
#' @param U: initialization of U or estimate from previous iteration
#' @param p: the weight vector for a component
#' @param S0: a covariance matrix in inverse-Wishart prior
#' @param n0: parameter in inverse-Wishart prior, n0 = \nu + R + 1
#' @param sigma2: initialization of the scalar value or estimate from previous iteration.
scaled_ed_reg_iid <- function(X, U, p, S0, n0, sigma2){
  m  <- ncol(X)
  I  <- diag(m)
  A  <- solve(U + I,U)
  B <- U %*% (I-A)
  bmat <- X %*% A
  U <- (crossprod(sqrt(p)*bmat)+sum(p)*B+sigma2*n0*S0)/(sum(p)+ n0)
  sigma2 <- m/sum(diag(S0 %*% solve(U)))
  return(list(U = U, sigma2 = sigma2))
}


#' Function to compute log posterior in scaled-ED model.
#' @param X: data matrix of size $n$ by $R$.
#' @param w: mixture component weights
#' @param U: prior covariance matrix
#' @param V: residual covariance matrix
#' @param S0: a covariance matrix in inverse-Wishart prior
#' @param n0: parameter in inverse-Wishart prior, n0 = \nu + R + 1
#' @param s: a vector of scalars for each prior covariance U.
compute_log_posterior_scaled_ed <- function(X, w, U, V, S0, n0, s){
  log_prior <- 0
  K <- length(w)
  for (k in 1:K){
    log_prior = log_prior + ldiwishart(W = U[,,k]/s[k], n0, S0)
  }
  loglik <- loglik_ud_iid_helper(X, w, U, V) 
  log_posterior <- loglik + log_prior
  return(log_posterior)
}

# Compute the log-density of the inverse Wishart at W with n - d - 1
# degrees of freedom and scale matrix n*S, ignoring terms that do not
# depend on X.
ldiwishart <- function (W, n, S)
  -n/2*(ldet(W) + tr(S %*% solve(W)))

# Compute the log-determinant of matrix A.
ldet <- function (A)
  as.numeric(determinant(A,logarithm = TRUE)$modulus)

# Compute trace of matrix A.
tr <- function (A)
  sum(diag(A))

#' Function to fit regularized ED. 
#' @param U: a 3d array contains k prior covariance matrix
em_fit_regularized_ed_scaled <- function(X, w, U, V, S0, n0, maxiter, tol, tol.obj){
  # initialize progress to store progress at each iteration
  
  progress <- as.data.frame(matrix(0, maxiter,5))
  names(progress) <- c("iter","loglik", "log_posterior", "delta.w","delta.u")
  progress$iter <- 1:maxiter

  # Get the number of samples (n) and the number of mixture components (k)
  n <- nrow(X)
  k <- length(w)
  m <- ncol(X)
  s <- rep(1, k)
  log_posterior0 = -Inf
  for (iter in 1:maxiter){
    # param values at previous iteration
    w0  <- w
    U0  <- U
    # E-step: calculate posterior probabilities using the current mu and sigmas
    P <- compute_posterior_probs_iid(X, w, U, V)
    # M-step: 
    # update covariance matrix 
    for (j in 1:k){
      res <- scaled_ed_reg_iid(X, U[,,j], P[,j], S0, n0, s[j])
      U[,,j] <- res$U
      s[j] <- res$sigma2
    }
    # update mixture weight
    w = colSums(P)/n
    # Compute log-likelihood.
    loglik = loglik_ud_iid_helper(X, w, U, V)
    log_posterior = compute_log_posterior_scaled_ed(X, w, U, V, S0, n0, s)
    dlogpost <- log_posterior - log_posterior0
    dw <- max(abs(w - w0))
    dU <- max(abs(U - U0))
    progress[iter,"loglik"]  <- loglik
    progress[iter, "log_posterior"] <- log_posterior
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    log_posterior0 <- log_posterior
    dparam  <- max(dw,dU)
    if (dparam < tol | dlogpost < log(1 + tol.obj))
      break
    }
    return(list(w = w, U = U, progress = progress[1:iter, ]))
}