#' @title Simulate Data from Ultimate Deconvolution Model
#'
#' @description Simulate data points from the Ultimate Deconvolution
#'   model. See \code{\link{ud_fit}} for the model definition.
#'
#' @param n Number of data points to simulate.
#' 
#' @param w A numeric vector of length k specifying the prior mixture
#'   weights. All entries must be non-negative, but need not sum to 1;
#'   the mixture weights are automatically normalized to sum to 1. If
#'   not provided, all the mixture weights are set to the same value. 
#'
#' @param U A list of length k specifying the covariance matrices in
#'   the mixture prior; list element \code{U[[i]]} is the m x m
#'   covariance matrix for the ith mixture component. For k = 1, U may
#'   also be a matrix.
#' 
#' @param V The m x m residual covariance matrix. If missing, \code{V}
#'   is set to the identity matrix.
#'
#' @return An n x m matrix in which each row is a draw from the
#'   Ultimate Deconvolution model.
#' 
#' @seealso \code{\link{ud_fit}}
#'
#' @importFrom mvtnorm rmvnorm
#' 
#' @export
#' 
simulate_ud_data <- function (n, w, U, V) {
      
  # Get the dimension of the data points (m) and the number of
  # mixture components (k).
  if (is.matrix(U))
    U <- list(U = U)
  m <- nrow(U[[1]])
  k <- length(U)
  if (n < 2)
    stop("n should be 2 or greater")
  if (m < 2)
    stop("The univariate case (m = 1) is not implemented")
  
  # Check the residual covariance matrix.
  if (missing(V))
    V <- diag(nrow(U[[1]]))
  if (!issemidef(V))
    stop("Input argument \"V\" should be a positive semi-definite matrix")
  
  # Check the prior covariance matrices.
  for (i in 1:k)
    if (!issemidef(U[[i]]))
      stop("All \"U\" matrices should be positive semi-definite")

  # Check the mixture weights.
  if (missing(w))
    w <- rep(1/k,k)
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\" should be a vector of length \"k\" ",
         "containing non-negative weights")
  w <- w/sum(w)

  # Draw mixture components according to the mixture weights.
  z <- sample(k,n,replace = TRUE,prob = w)
  
  # Draw the data points.
  X <- matrix(0,n,m)
  for (j in 1:k) {
    i <- which(z == j)
    if (length(i) > 0) {
      T <- V + U[[j]]
      X[i,] <- rmvnorm(length(i),sigma = T)
    }
  }
  rownames(X) <- paste0("s",1:n)
  colnames(X) <- rownames(V)
  return(X)
}


# Function to simulate data from EBNM model for one component.
# It first simulate the true mean and then observed data. 
# @param n Number of data points to simulate.
# @param U Prior covariance matrix
# @param V Residual covariance matrix
simulate_ebmn_data <- function(n, U, V){
  m = ncol(U)
  X = matrix(NA, nrow = n, ncol = m)
  theta <- rmvnorm(n, sigma = U)
  for (i in 1:n)
    X[i, ] <- rmvnorm(1, theta[i, ], sigma = V)
  return(list(X = X, theta = theta))
}

