ed_reg_iid <- function(X, U, p){
  m  <- ncol(X)
  I  <- diag(m)
  A  <- solve(U + I,U)
  B <- U %*% (I-A)
  bmat <- X %*% A
  Phi <- 2*m*diag(m)
  U <- (crossprod(sqrt(p)*bmat)+sum(p)*B+Phi)/(sum(p)+ 2*m)
}


#' Function to fit regularized ED. 
#' @param U: a 3d array contains k prior covariance matrix
em_fit_regularized_ed <- function(X, w, U, V, maxiter, tol, tol.lik){
  # initialize progress to store progress at each iteration
  
  progress <- as.data.frame(matrix(0, maxiter,4))
  names(progress) <- c("iter","loglik","delta.w","delta.u")
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
      U[,,j] = ed_reg_iid(X, U[,,j], P[,j])
    }
    # update mixture weight
    w = colSums(P)/n
    # Compute log-likelihood.
    loglik = loglik_ud_iid_helper(X, w, U, V)
  
    dw <- max(abs(w - w0))
    dU <- max(abs(U - U0))
    dloglik <- loglik - loglik0
    progress[iter,"loglik"]  <- loglik
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    
    dparam  <- max(dw,dU)
    loglik0 <- loglik
    if (dparam < tol | dloglik < log(1 + tol.lik))
      break
    }
    return(list(w = w, U = U, progress = progress[1:iter, ]))
}