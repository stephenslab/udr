# A test that mvebnm also works for the univariate case.

# SIMULATE DATA
# -------------
n <- 800
S <- 2
U <- list(0,0.1,1,10)
w <- c(0.8,0.1,0.075,0.025)

# Add row and column names to the variables.
k           <- length(w)
names(w)    <- paste0("k",1:k)
names(U)    <- paste0("k",1:k)
  
# Simulate draws from the normal means model.
X        <- simulate_mvebnm_data(n,w,U,S)
names(X) <- paste0("s",1:n)

# FIT MODEL
# ---------
set.seed(1)
numiter <- 100
t1 <- system.time(
  fit1 <- mvebnm(X,k = k,S = S,verbose = TRUE,
                 control = list(version = "R",update.U = "ed",
                                update.S = "em",maxiter = numiter)))

set.seed(1)
t2 <- system.time(
  fit2 <- mvebnm(X,k = k,S = S,verbose = TRUE,
                 control = list(version = "Rcpp",update.U = "ed",
                                update.S = "em",maxiter = numiter)))

print(fit1$loglik - fit2$loglik)
print(range(fit1$progress$loglik - fit2$progress$loglik))
print(range(fit1$w - fit2$w))
print(fit1$S - fit2$S)
print(range(unlist(fit1$U) - unlist(fit2$U)))

set.seed(1)
t3 <- system.time(
  fit3 <- mvebnm(X,k = k,S = S,verbose = TRUE,
                 control = list(version = "R",maxiter = numiter,tol = 0,
                                update.S = "em")))
set.seed(1)
t4 <- system.time({
  fit4 <- mvebnm(X,k = k,S = S,verbose = TRUE,
                 control = list(version = "Rcpp",tol = 0,maxiter = 50,
                                update.S = "em"))
  fit4 <- mvebnm(X,fit0 = fit4,verbose = TRUE,
                 control = list(version = "Rcpp",tol = 0,maxiter = 50,
                                update.S = "em"))
})

print(fit3$loglik - fit4$loglik)
print(range(fit3$progress$loglik - fit4$progress$loglik))
print(range(fit3$w - fit4$w))
print(fit3$S - fit4$S)
print(range(unlist(fit3$U) - unlist(fit4$U)))

y1 <- fit1$progress$loglik
y2 <- fit2$progress$loglik
y3 <- fit3$progress$loglik
y4 <- fit4$progress$loglik
y  <- max(c(y1,y2,y3,y4))
plot(1:numiter,y - y1 + 0.01,col = "dodgerblue",type = "l",log = "y",lwd = 2,
     xlab = "iteration",ylab = "dist. from best loglik",ylim = c(0.01,500))
lines(1:numiter,y - y2 + 0.01,col = "darkblue",lwd = 2,lty = "dashed")
lines(1:numiter,y - y3 + 0.01,col = "magenta",lwd = 2)
lines(1:numiter,y - y4 + 0.01,col = "gold",lwd = 2,lty = "dashed")
