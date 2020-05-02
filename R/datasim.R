# Simulate n data points from the "extreme deconvolution" mixture
# model: each data point is drawn from a mixture of normals in which
# each normal mixture component k has zero mean and covariance matrix
# S + U[[k]]. Input argument w specifies the mixture weights.
#'
#' @export
#' 
datasim.ed <- function (n, w, U, S) {

  # Get the number of outcomes (m) and the number of mixture
  # components (k).
  m <- nrow(S)
  k <- length(w)

  # Check that all the covariance matrices are s.p.d.
  for (j in 1:k)
    tryCatch(R <- chol(S + U[[j]]),
             error = function (e)
               stop(paste("One or more covariance matrices S + U are not",
                          "symmetric positive definite")))
  
  # Initialize the n x m output matrix.
  X <- matrix(0,n,m)

  # Repeat for each sample.
  for (i in 1:n) {

    # Draw a mixture component according to the mixture weights.
    j <- sample(k,1,prob = w)
    
    # Draw the data point.
    T     <- S + U[[j]]
    X[i,] <- rmvnorm(1,sigma = T)
  }

  return(X)
}
