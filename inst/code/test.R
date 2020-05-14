# TO DO: Explain here what this script does, and how to use it.
library(mvtnorm)

# SIMULATE DATA
# -------------
# Simulate data points x ~ N(z,S), in which z ~ N(0,V).
set.seed(1)
n <- 1000
S <- rbind(c(0.8,0.4),
           c(0.4,1.2))
V <- rbind(c(1.0,0.9),
           c(0.9,1.0))
Z <- rmvnorm(n,sigma = V)
X <- matrix(0,n,2)
for (i in 1:n)
  X[i,] <- rmvnorm(1,Z[i,],S)

# ESTIMATE COVARIANCE, S
# ----------------------
S      <- diag(2)
loglik <- rep(0,40)
for (iter in 1:40) {

  # Update S.
  Snew <- matrix(0,2,2)
  for (i in 1:n) {
    x    <- X[i,]
    S1   <- solve(V %*% solve(S) + diag(2)) %*% V
    mu1  <- drop(S1 %*% solve(S,x))
    Snew <- Snew + S1 + tcrossprod(x - mu1)
  }
  Snew <- Snew/n

  # Compute the log-likelihood at the new estimate of S.
  S <- Snew
  loglik[iter] <- sum(dmvnorm(X,sigma = S + V,log = TRUE))
}
