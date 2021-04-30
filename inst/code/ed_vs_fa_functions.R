# TO DO: Explain here what this function does, and how to use it.
ed <- function (X, U, V, numiter = 100) {
  n <- nrow(X)
  m <- ncol(X)
  
  # These two variables are used to keep track of the algorithm's
  # progress.
  loglik  <- rep(0,numiter)
  maxdiff <- rep(0,numiter)
  
  # Iterate the EM updates.
  for (iter in 1:numiter) {

    # Save the current parameter estimates.
    U0 <- U
    V0 <- V

    # E STEP
    # ------
    # Compute the posterior means and covariances of the latent v's.
    S  <- solve(solve(U) + solve(V))
    mu <- t(S %*% solve(V,t(X)))
    
    # M STEP
    # ------
    # Update the prior covariance matrix, U.
    U <- S + crossprod(mu)/n
    
    # Update the residual covariance matrix, V.
    V <- S + crossprod(X - mu)/n
    
    # Record the algorithm's progress.
    T <- U + V
    loglik[iter]  <- sum(dmvnorm(X,sigma = T,log = TRUE))
    maxdiff[iter] <- max(max(abs(U - U0)),max(abs(V - V0)))
  }

  return(list(U = U,V = V,loglik = loglik,maxdiff = maxdiff))
}

# TO DO: Explain here what this function does, and how to use it.
fa <- function (X, U, V, numiter = 100) {

}
