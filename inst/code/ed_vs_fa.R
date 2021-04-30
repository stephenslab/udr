# TO DO: Explain here what this script is for, and how to use it.
library(mvtnorm)
source("ed_vs_fa_functions.R")

# Simulate data.
set.seed(1)
n <- 200
U <- rbind(c(1.0,0.9),
           c(0.9,1.0))
V <- rbind(c(0.8,0.2),
           c(0.2,1.5))
T <- U + V
X <- rmvnorm(n,sigma = T)

# Perform the EM updates derived using the ED data augmentation.
numiter <- 20
fit1 <- ed(X,U,V,numiter = numiter)

# Perform the EM updates derived using the FA data augmentation.
# TO DO

# Plot the improvement in the solution over time.
y1 <- fit1$loglik
y1 <- max(y1) - y1 + 1e-10
plot(1:numiter,y1,type = "l",lwd = 2,col = "darkorange",log = "y",
     xlab = "iteration",ylab = "loglik diff")
