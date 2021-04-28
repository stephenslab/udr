set.seed(1)
n <- 400
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
rownames(V) <- c("d1","d2")
colnames(V) <- c("d1","d2")
X <- simulate_ud_data(n,w,U,V)
fit0 <- ud_init(X,
                U_scaled = U,
                U_rank1 = list(tcrossprod(c(-1,1)),tcrossprod(c(1,2))),
                n_unconstrained = 2,
                V = V) # V = rep(list(V),n)
control <- list(maxiter = 100,
                resid.update = "em",
                scaled.update = "none",
                rank1.update = "none",
                unconstrained.update = "ed",
                version = "Rcpp")
fit1 <- ud_fit(fit0,control = control)
y <- fit1$progress$loglik
y <- max(y) - y + 0.01
plot(1:100,y,col = "darkblue",type = "l",lwd = 2,log = "y",
     xlab = "iteration",ylab = "loglik difference")

