# TO DO: Explain here what this script is for, and how to use it.
library(mvtnorm)
source("ed_vs_fa_functions.R")

# Simulate data.
set.seed(1)
n <- 200
U <- rbind(c(1.0,0.5),
           c(0.5,1.0))
V <- rbind(c(0.8,0.2),
           c(0.2,1.5))
T <- U + V
X <- rmvnorm(n,sigma = T)

# Perform the EM updates derived using the ED data augmentation.
numiter <- 40
fit1 <- ed(X,U,V,numiter)

# Perform the EM updates derived using the FA data augmentation.
fit2 <- fa(X,U,V,numiter)

# Plot the improvement in the solution over time.
e  <- 1e-10
y1 <- fit1$loglik
y2 <- fit2$loglik
y  <- c(y1,y2)
y1 <- max(y) - y1 + e
y2 <- max(y) - y2 + e
plot(1:numiter,y1,type = "l",lwd = 2,col = "darkorange",log = "y",
     xlab = "iteration",ylab = "loglik diff",ylim = c(e,max(c(y1,y2))))
lines(1:numiter,y2,col = "dodgerblue",lwd = 2,lty = "dashed")
