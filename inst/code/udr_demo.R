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

# FIT MODEL
# ---------
# This is the simplest invocation of ud_init and ud_fit.
set.seed(1)
fit1 <- ud_init(X)
fit1 <- ud_fit(fit1,control = list(maxiter = 10))
print(summary(fit1))
plot(fit1$progress$iter,max(fit1$progress$loglik) - fit1$progress$loglik,
     type = "l",col = "dodgerblue",lwd = 2,log = "",xlab = "iteration",
     ylab = "dist to best loglik")

# This is a more complicated invocation of ud_init that overrides some
# of the defaults.
fit2 <- ud_init(X,U_scaled = U,n_rank1 = 1,n_unconstrained = 1,V = V)
fit2 <- ud_fit(fit2,control = list(maxiter = 10))
print(summary(fit2))
plot(fit2$progress$iter,max(fit2$progress$loglik) - fit2$progress$loglik,
     type = "l",col = "dodgerblue",lwd = 2,log = "",xlab = "iteration",
     ylab = "dist to best loglik")
