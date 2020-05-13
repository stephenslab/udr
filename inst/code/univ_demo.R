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
