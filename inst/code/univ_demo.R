# A test that mvebnm also works for the univariate case.

# SIMULATE DATA
# -------------
n <- 400
S <- 2
U <- list(0,0.1,1,10)
w <- c(0.8,0.1,0.075,0.025)

# Add row and column names to the variables.
k           <- length(w)
names(w)    <- paste0("k",1:k)
names(U)    <- paste0("k",1:k)
  
# Simulate draws from the normal means model.
X        <- simulate_ebmvnm_data(n,w,U,S)
names(X) <- paste0("s",1:n)

# FIT MODEL
# ---------
# TO DO: Explain these lines of code in greater detail.
set.seed(1)
numiter <- 100
t1 <- system.time(
  fit1 <- mvebnm(X,k = k,S = S,control = list(version = "R",update.U = "em",
                                   maxiter = numiter)))

set.seed(1)
t2 <- system.time(
  fit2 <- mvebnm(X,k = k,S = S,control = list(version = "Rcpp",update.U = "em",
                                              maxiter = numiter)))

print(fit1$loglik - fit2$loglik)
print(range(fit1$progress$loglik - fit2$progress$loglik))
print(range(fit1$w - fit2$w))
print(fit1$S - fit2$S)
print(range(unlist(fit1$U) - unlist(fit2$U)))

set.seed(1)
t3 <- system.time(
  fit3 <- mvebnm(X,k = k,S = S,control = list(version = "R",
                                              maxiter = numiter)))
set.seed(1)
t4 <- system.time(
  fit4 <- mvebnm(X,k = k,S = S,control = list(version = "Rcpp",

                                              maxiter = numiter)))
print(fit3$loglik - fit4$loglik)
print(range(fit3$progress$loglik - fit4$progress$loglik))
print(range(fit3$w - fit4$w))
print(fit3$S - fit4$S)
print(range(unlist(fit3$U) - unlist(fit4$U)))
