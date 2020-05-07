# This is a small demo script to verify initial implementation of the
# mvebnm methods.
set.seed(1)

# TO DO: Explain what these next lines of code do.
n <- 800
S <- rbind(c(0.8,0.4),
           c(0.4,1.2))
U <- list(none   = rbind(c(1e-6,0),
                         c(0,1e-6)),
          shared = rbind(c(1.0,0.9),
                         c(0.9,1.0)),
          only1  = rbind(c(1,0),
                         c(0,1e-6)),
          only2  = rbind(c(1e-6,0),
                         c(0,1)))
w <- c(0.8,0.1,0.075,0.025)

# TO DO: Explain what these next lines of code do.
rownames(S) <- c("d1","d2")
colnames(S) <- c("d1","d2")
k <- length(U)
for (i in 1:k) {

}
  
# TO DO: Explain what these next lines of code do.
X <- simulate_ebmvnm_data(n,w,U,S)
rownames(X) <- paste0("s",1:n)
