# Estimate the matrices U and V in the model x ~ N(0,U + V) using the
# EM algorithm derived using the Extreme Deconvolution (ED) data
# augmentation x ~ N(v,V), v ~ N(0,U).
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
    # Compute the posterior means (mu) and covariances (S) of the
    # latent variables v.
    S  <- solve(solve(U) + solve(V))
    mu <- t(S %*% solve(V,t(X)))
    
    # M STEP
    # ------
    # Update the residual covariance matrix, V.
    V <- S + crossprod(X - mu)/n
    
    # Update the prior covariance matrix, U.
    U <- S + crossprod(mu)/n
    
    # Record the algorithm's progress.
    T <- U + V
    loglik[iter]  <- sum(dmvnorm(X,sigma = T,log = TRUE))
    maxdiff[iter] <- max(max(abs(U - U0)),max(abs(V - V0)))
  }

  return(list(U = U,V = V,loglik = loglik,maxdiff = maxdiff))
}

# Estimate the matrices U and V in the model x ~ N(0,U + V) using the
# EM algorithm derived using the factor analysis (FA) data
# augmentation x ~ N(L*z,V), z ~ N(0,I) such that tcrossprod(L) = U.
fa <- function (X, U, V, numiter = 100) {
  n <- nrow(X)
  m <- ncol(X)
  L <- t(chol(U))

  # These two variables are used to keep track of the algorithm's
  # progress.
  loglik  <- rep(0,numiter)
  maxdiff <- rep(0,numiter)
  
  # Iterate the EM updates.
  for (iter in 1:numiter) {

    # Save the current parameter estimates.
    L0 <- L
    V0 <- V

    # E STEP
    # ------
    # Compute the posterior means (mu) and covariances (S) of the
    # latent variables z. Here, L is the left Cholesky factor of U so
    # that that tcrossprod(L) = U.
    I  <- diag(m)
    S  <- I - t(L) %*% solve(U + V,L)
    mu <- t(solve(U + V,t(X))) %*% L

    # Compute the posterior means (mv) and covariances (Sv) of t he
    # latent variables v = L*z.
    mv <- mu %*% t(L)
    Sv <- L %*% S %*% t(L)
    
    # M STEP
    # ------
    # Update V.
    V <- Sv + crossprod(X - mv)/n

    # Update L and U.
    L <- t(solve(n*S + crossprod(mu),t(crossprod(X,mu))))
    U <- tcrossprod(L)
    
    # Record the algorithm's progress.
    T <- U + V
    loglik[iter]  <- sum(dmvnorm(X,sigma = T,log = TRUE))
    maxdiff[iter] <- max(max(abs(L - L0)),max(abs(V - V0)))
  }
  
  return(list(U = U,L = L,V = V,loglik = loglik,maxdiff = maxdiff))
}

