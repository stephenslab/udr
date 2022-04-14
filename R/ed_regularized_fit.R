
#' Function for regularized ED where we specify an inverse Wishart
#' prior on U. U \sim W_R^{-1}(\Phi, \nu). Default parameter values:
#' \nu = R and \Phi = (\nu+R+1)*I.
#' The derivation for U is available in overleaf eq.(33).
#' https://www.overleaf.com/read/vrgwpskkhbpj
ed_reg_iid <- function(X, U, p, Phi = (nu+ncol(U)+1)*diag(ncol(U)), nu = diag(ncol(U))){
  m  <- ncol(X)
  I  <- diag(m)
  A  <- solve(U + I,U)
  B <- U %*% (I-A)
  bmat <- X %*% A
  U <- (crossprod(sqrt(p)*bmat)+sum(p)*B+Phi)/(sum(p)+ nu + m + 1)
}


#' Function to compute log posterior. 
compute_log_posterior_ed <- function(X, w, U, V, Phi, nu){
  log_prior <- 0
  K <- length(w)
  for (k in 1:K){
    log_prior = log_prior + dinvwishart((U[,,k]+t(U[,,k]))/2, nu = nu, S = Phi, log=TRUE)
  }
  loglik <- loglik_ud_iid_helper(X, w, U, V) 
  log_posterior <- loglik + log_prior
  return(log_posterior)
}


#' Function to fit regularized ED. 
#' @param U: a 3d array contains k prior covariance matrix
em_fit_regularized_ed <- function(X, w, U, V, Phi, nu, maxiter, tol, tol.lik){
  # initialize progress to store progress at each iteration
  
  progress <- as.data.frame(matrix(0, maxiter,5))
  names(progress) <- c("iter","loglik", "log_posterior", "delta.w","delta.u")
  progress$iter <- 1:maxiter

  # Get the number of samples (n) and the number of mixture components (k)
  n <- nrow(X)
  k <- length(w)
  m <- ncol(X)
  for (iter in 1:maxiter){
    # param values at previous iteration
    w0  <- w
    U0  <- U
    loglik0 = -Inf
    # E-step: calculate posterior probabilities using the current mu and sigmas
    P <- compute_posterior_probs_iid(X, w, U, V)
    # M-step: 
    # update covariance matrix 
    for (j in 1:k){
      U[,,j] = ed_reg_iid(X, U[,,j], P[,j], Phi, nu)
    }
    # update mixture weight
    w = colSums(P)/n
    # Compute log-likelihood.
    loglik = loglik_ud_iid_helper(X, w, U, V)
    log_posterior = compute_log_posterior(X, w, U, V, Phi, nu)

    dw <- max(abs(w - w0))
    dU <- max(abs(U - U0))
    dloglik <- loglik - loglik0
    progress[iter,"loglik"]  <- loglik
    progress[iter, "log_posterior"] <- log_posterior
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    
    dparam  <- max(dw,dU)
    loglik0 <- loglik
    if (dparam < tol | dloglik < log(1 + tol.lik))
      break
    }
    return(list(w = w, U = U, progress = progress[1:iter, ]))
}