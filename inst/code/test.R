# TO DO: Explain here what this script does, and how to use it.
library(mvtnorm)

# SIMULATE DATA
# -------------
set.seed(1)
n <- 100
S <- rbind(c(0.8,0.4),
           c(0.4,1.2))
V <- rbind(c(1.0,0.9),
           c(0.9,1.0))
Z <- rmvnorm(n,sigma = V)
X <- matrix(0,n,2)
for (i in 1:n)
  X[i,] <- rmvnorm(1,Z[i,],S)
