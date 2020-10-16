# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function (x)
  expect_equal(diff(x) >= 0,rep(TRUE,length(x) - 1))

# Simulate data from the (univariate) normal means model.
simulate_ebnm_data <- function (n) {

  # These parameters specify the model used to simulate the data.
  S <- 2
  U <- list(0,0.1,1,10)
  w <- c(0.8,0.1,0.075,0.025)

  # Add row and column names to the variables.
  k        <- length(w)
  names(w) <- paste0("k",1:k)
  names(U) <- paste0("k",1:k)
  
  # Simulate draws from the (univariate) normal means model.
  X        <- simulate_ud_data(n,w,U,S)
  names(X) <- paste0("s",1:n)

  # Output the data, and the parameters of the model used to simulate
  # the data.
  return(list(X = X,w = w,U = U,S = S))
}

# Simulate bivariate data points drawn from a mixture of four normals.
simulate_ud_data_2d <- function (n) {
    
  # These parameters specify the model used to simulate the data.
  S <- rbind(c(0.8,0.2),
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
  rownames(S) <- c("d1","d2")
  colnames(S) <- c("d1","d2")
  k           <- length(w)
  names(w)    <- names(U)
  for (i in 1:k) {
    rownames(U[[i]]) <- c("d1","d2")
    colnames(U[[i]]) <- c("d1","d2")
  }

  # Simulate draws from the multivariate normal means model.
  X           <- simulate_ud_data(n,w,U,S)
  rownames(X) <- paste0("s",1:n)

  # Output the data, and the parameters of the model used to simulate
  # the data.
  return(list(X = X,w = w,U = U,S = S))
}
