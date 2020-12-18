#' @title Simulate Data from Ultimate Deconvolution Model
#'
#' @description Simulate data points from the Ultimate Deconvolution
#'   model. The univariate case (m = 1) is also handled. See
#'   \code{\link{ud_fit}} for the model definition.
#'
#' @param n Number of data points to simulate.
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
#' @param V The m x m residual covariance matrix.
#'
#' @return An n x m matrix in which each row is a draw from the
#'   Ultimate Deconvolution model. For the univariate case (m = 1), a
#'   vector is returned.
#' 
#' @seealso \code{\link{ud_fit}}
#'
#' @importFrom mvtnorm rmvnorm
#' 
#' @export
#' 
simulate_ud_data <- function (n, w, U, V) {
      
  # Check the residual covariance matrix, V.
  V <- as.matrix(V)
  if (!issemidef(V))
    stop("Input argument \"V\" should be a positive semi-definite matrix")
  
  # Check the prior covariance matrices, U.
  U <- lapply(U,as.matrix)
  if (!(is.list(U) && verify_prior_covariances(U,V)))
    stop("Input argument \"U\" should be list in which each list element ",
         "U[[i]] is a (symmetric) positive semi-definite matrix, and",
         "V + U[[i]] is symmetric positive definite")

  # Get the dimension of the data points (m) and the number of
  # mixture components (k).
  m <- nrow(V)
  k <- length(U)

  # Check the mixture weights, w.
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\" should be a vector of length \"k\" ",
         "containing non-negative weights")
  w <- w/sum(w)

  # Initialize storage for the data points.
  X <- matrix(0,n,m)

  # Draw mixture components according to the mixture weights.
  z <- sample(k,n,replace = TRUE,prob = w)
  
  # Draw the data points.
  for (j in 1:k) {
    i     <- which(z == j)
    T     <- V + U[[j]]
    X[i,] <- rmvnorm(length(i),sigma = T)
  }
  colnames(X) <- rownames(V)
  return(drop(X))
}

