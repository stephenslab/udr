#' @title Simulate Data from Multivariate Normal Means Model
#'
#' @description Simulate n data points from the multivariate normal
#'   means model. See \code{\link{mvebnm}} for the model definition.
#'
#' @param w A numeric vector of length k specifying the prior mixture
#'   weights. All entries must be non-negative, but need not sum to 1;
#'   the mixture weights are automatically normalized to sum to 1.
#'
#' @param U A list of length k specifying the covariance matrices in
#'   the mixture-of-multivariate-normals prior; list element
#'   \code{U[[i]]} is the m x m covariance matrix for the ith mixture
#'   component.
#' 
#' @param S The m x m residual covariance matrix.
#'
#' @return Describe output here.
#' 
#' @seealso \code{\link{mvebnm}}
#' 
#' @export
#' 
simulate_ebmvnm_data <- function (n, w, U, S) {

  # Get the number of outcomes (m) and the number of mixture
  # components (k).
  m <- nrow(S)
  k <- length(w)

  # Check that all the (marginal) covariance matrices are s.p.d.
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
