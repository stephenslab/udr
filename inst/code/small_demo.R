# This is a small demo script to verify initial implementation of the
# mvebnm methods.
set.seed(1)

# SIMULATE DATA
# -------------
# These variables specify how the data are simulated: n is the number
# of samples drawn from the multivariate normal means model,
#
#   w[1]*N(0,S + U[[1]]) + ... + w[4]*N(0,S + U[[4]]).
#
n <- 2000
S <- rbind(c(0.8,0.4),
           c(0.4,1.2))
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
rownames(S) <- c("d1","d2")
colnames(S) <- c("d1","d2")
k           <- length(w)
names(w)    <- paste0("k",1:k)
for (i in 1:k) {
  rownames(U[[i]]) <- c("d1","d2")
  colnames(U[[i]]) <- c("d1","d2")
}
  
# Simulate draws from the multivariate normal means model.
X           <- simulate_ebmvnm_data(n,w,U,S)
rownames(X) <- paste0("s",1:n)

# FIT MODEL
# ---------
# TO DO: Explain these lines of code in greater detail.
set.seed(1)
fit1 <- mvebnm(X,k = k,S = S,control = list(version = "R",maxiter = 40))
y <- fit1$progress$loglik
plot(fit1$progress$iter,max(y) - y + 0.01,col = "dodgerblue",type = "l",
     log = "y",lwd = 2,xlab = "iteration",ylab = "dist. from best loglik")
# fit1 <- mvebnm(X,k,control = list(update.U = "em",version = "R"))
# cat("\n")
# set.seed(1)
# fit2 <- mvebnm(X,k,control = list(update.U = "em",version = "Rcpp"))
# print(fit1$loglik - fit2$loglik)
