# TO DO: Explain here what this function does, and how to use it.
ed <- function (X, U, V, numiter = 100) {
  n <- nrow(X)
  m <- ncol(X)
  
  # These two variables are used to keep track of the algorithm's
  # progress.
  loglik <- rep(0,numiter)
  maxd   <- rep(0,numiter)
  
  # Iterate the EM updates.
  for (iter in 1:numiter) {

    # Save the current parameter estimates.
    U0 <- U
    V0 <- V

    # E STEP
    # ------
    mu <- matrix(0,n,m)
    S  <- solve(solve(U) + solve(V))
    for (i in 1:n)
      mu[i,] <- S %*% solve(V,X[i,])

    # M STEP
    # ------
    # Update the prior covariance matrix, U.
    # TO DO.

    # Update the residual covariance matrix, V.
    V <- n*S
    for (i in 1:n)
      V <- V + tcrossprod(X[i,] - mu[i,])
    V <- V/n

    # Record the algorithm's progress.
    T <- U + V
    loglik[iter] <- sum(dmvnorm(X,sigma = T,log = TRUE))
    maxd[iter] <- max(max(abs(U - U0)),max(abs(V - V0)))
  }

  return(list(U = U,V = V,loglik = loglik,maxd = maxd))
}

# TO DO: Explain here what this function does, and how to use it.
fa <- function (X, U, V, numiter = 100) {

}
