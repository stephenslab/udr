# Small script to motivate the parameterization U = V*U'*V.
library(mvtnorm)
set.seed(1)

# Try with different penalty strengths.
a <- 0.01

# Simulate 2-d data points, x ~ N(0,S).
n <- 200
U <- diag(2)
V <- diag(c(10,1))
S <- V %*% U %*% V
X <- rmvnorm(n,sigma = S)

# Get the maximum-likelihood estimate of S.
Smle <- cov(X)

# Encourage the eigenvalues of S to be similar ("well-conditioned").
getcov1 <- function (x) {
  S      <- matrix(0,2,2)
  x[1:2] <- abs(x[1:2])
  S[1:4] <- x[c(1,3,3,2)]
  return(S)
}

f1 <- function (vars) {
  S  <- getcov1(vars)
  e  <- eigen(S)$values
  ll <- dmvnorm(X,sigma = S,log = TRUE)
  return(-mean(ll) + a*abs(diff(e)))
}

res <- optim(c(1,1,0.5),f1,method = "Nelder-Mead")
S1  <- getcov1(res$par)

# Encourage the eigenvalues of U to be similar, where S = V*U*V.
getcov2 <- function (x) {
  x[1:4] <- abs(x[1:4])
  V      <- diag(x[1:2])
  U      <- matrix(0,2,2)
  U[1:4] <- x[c(3,5,5,4)]
  S      <- V %*% U %*% V
  return(list(U = U,V = V,S = S))
}
    
f2 <- function (vars) {
  res <- getcov2(vars)
  ll  <- dmvnorm(X,sigma = res$S,log = TRUE)
  e   <- eigen(res$U)$values
  return(-mean(ll) + a*abs(diff(e)))
}
    
res <- optim(c(rep(1,4),0.5),f2,method = "Nelder-Mead")
S2  <- getcov2(res$par)$S

# Compare the estimates with the ground-truth covariance (S).
cat("S =\n")
print(round(S,digits = 3))
cat("Smle =\n")
print(round(Smle,digits = 3))
cat("S1 =\n")
print(round(S1,digits = 3))
cat("S2 =\n")
print(round(S2,digits = 3))
