# This is a small demo script to verify initial implementation of the
# Ultimate Deconvolution methods.
set.seed(1)

# SIMULATE DATA
# -------------
# These variables specify how the data are simulated: n is the number
# of samples drawn from the multivariate normal means model,
#
#   w[1]*N(0,V + U[[1]]) + ... + w[4]*N(0,V + U[[4]]).
#
n <- 4000
V <- rbind(c(0.8,0.2),
           c(0.2,1.5))
U <- list(none   = rbind(c(0,0),
                         c(0,0)),
          shared = rbind(c(1.0,0.9),
                         c(0.9,1.0)),
          only1  = rbind(c(1,0),
                         c(0,0)),
          only2  = rbind(c(0,0),
                         c(0,1)))
w <- c(0.8,0.1,0.075,0.025)

# Add row and column names to the variables.
rownames(V) <- c("d1","d2")
colnames(V) <- c("d1","d2")
k           <- length(w)
names(w)    <- paste0("k",1:k)
for (i in 1:k) {
  rownames(U[[i]]) <- c("d1","d2")
  colnames(U[[i]]) <- c("d1","d2")
}
  
# Simulate draws from the multivariate normal means model.
X           <- simulate_ud_data(n,w,U,V)
rownames(X) <- paste0("s",1:n)

# INITIALIZE FIT
# --------------
# This is the simplest invocation of ud_init.
set.seed(1)
fit1 <- ud_init(X)
print(summary(fit1))

# This is a more complicated invocation of ud_init that overrides some
# of the defaults.
fit2 <- ud_init(X,U_scaled = U,n_rank1 = 1,n_unconstrained = 1,V = V)
print(summary(fit2))

# FIT MODEL
# ---------
set.seed(1)
numiter <- 100
t1 <- system.time(
  fit1 <- ud_fit(X,k = k,S = S,
                 control = list(version = "R",update.U = "ed",
                                tol = 0,maxiter = numiter)))
set.seed(1)
t2 <- system.time(
  fit2 <- ud_fit(X,k = k,S = S,
                 control = list(version = "Rcpp",update.U = "ed",
                                tol = 0,maxiter = numiter)))
set.seed(1)
t3 <- system.time(
  fit3 <- ud_fit(X,k = k,S = S,
                 control = list(version = "R",update.U = "teem",
                                tol = 0,maxiter = numiter)))
print(lapply(fit3$U,function (x) eigen(x)$values))
set.seed(1)
t4 <- system.time(
  fit4 <- ud_fit(X,k = k,S = S,
                 control = list(version = "Rcpp",update.U = "teem",
                                tol = 0,maxiter = 50)))

fit4a <- ud_fit(X,fit0 = fit4,
               control = list(version = "Rcpp",update.U = "teem",
                              tol = 0,maxiter = 50))

fit4b <- ud_fit(X,fit0 = fit4,
               control = list(version = "Rcpp",update.U = "teem",
                              update.w = "mixsqp",tol = 0,
                              maxiter = 50))

set.seed(1)
fit5 <- ud_fit(X,k = k,S = S,
               control = list(version = "R",update.U = "teem",
                              update.S = "em",tol = 0,
                              maxiter = numiter))

set.seed(1)
fit6 <- ud_fit(X,k = k,S = S,
               control = list(version = "Rcpp",update.U = "teem",
                              update.S = "em",tol = 0,
                              maxiter = numiter))

print(t1)
print(t2)
print(t3)
print(t4)

print(fit1$loglik - fit2$loglik)
print(range(fit1$progress$loglik - fit2$progress$loglik))
print(range(fit1$w - fit2$w))
print(range(fit1$S - fit2$S))
for (i in 1:k)
  print(range(fit1$U[[i]] - fit2$U[[i]]))

print(fit3$loglik - fit4$loglik)
print(range(fit3$progress$loglik - fit4$progress$loglik))
print(range(fit3$w - fit4$w))
print(range(fit3$S - fit4$S))
for (i in 1:k)
  print(range(fit3$U[[i]] - fit4$U[[i]]))

y1  <- fit1$progress$loglik
y2  <- fit2$progress$loglik
y3  <- fit3$progress$loglik
y4a <- fit4a$progress$loglik
y4b <- fit4b$progress$loglik
y5  <- fit5$progress$loglik
y6  <- fit6$progress$loglik
y   <- max(c(y1,y2,y3,y4a,y4b,y5,y6))
plot(1:numiter,y - y1 + 0.01,col = "dodgerblue",type = "l",log = "y",lwd = 2,
     xlab = "iteration",ylab = "dist. from best loglik",ylim = c(0.01,500))
lines(1:numiter,y - y2 + 0.01,col = "darkblue",lwd = 2,lty = "dashed")
lines(1:numiter,y - y3 + 0.01,col = "magenta",lwd = 2)
lines(1:numiter,y - y4a + 0.01,col = "gold",lwd = 2,lty = "dashed")
lines(1:numiter,y - y4b + 0.01,col = "darkorange",lwd = 2,lty = "dashed")
lines(1:numiter,y - y5 + 0.01,col = "limegreen",lwd = 2)
lines(1:numiter,y - y5 + 0.01,col = "seagreen",lwd = 2,lty = "dashed")
