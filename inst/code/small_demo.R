# Simulate data from a UD model.
set.seed(1)
n <- 1000
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

# This is the simplest invocation of ud_init.
fit1 <- ud_init(X)
fit1 <- ud_fit(fit1)
summary(fit1)
plot(fit1$progress$iter,
     max(fit1$progress$loglik) - fit1$progress$loglik + 0.1,
     type = "l",col = "dodgerblue",lwd = 2,log = "y",xlab = "iteration",
     ylab = "dist to best loglik")
