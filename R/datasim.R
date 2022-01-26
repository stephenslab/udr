#' @title Simulate Data from Ultimate Deconvolution Model
#'
#' @description Simulate data points from the Ultimate Deconvolution
#'   model. See \code{\link{ud_fit}} for the model definition.
#'
#' @param n Number of data points to simulate.
#' 
#' @param w A numeric vector of length k specifying the prior mixture
#'   weights. All entries must be non-negative, but need not sum to 1;
#'   the mixture weights are automatically normalized to sum to 1.
#'
#' @param U A list of length k specifying the covariance matrices in
#'   the mixture prior; list element \code{U[[i]]} is the m x m
#'   covariance matrix for the ith mixture component.
#' 
#' @param V The m x m residual covariance matrix.
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
  m <- nrow(V)
  k <- length(U)
  if (n < 2)
    stop("n should be 2 or greater")
  if (m < 2)
    stop("The univariate case (m = 1) is not implemented")
  
  # Check the residual covariance matrix.
  if (!issemidef(V))
    stop("Input argument \"V\" should be a positive semi-definite matrix")
  
  # Check the prior covariance matrices.
  for (i in 1:k)
    if (!issemidef(U[[i]]))
      stop("All \"U\" matrices should be positive semi-definite")

  # Check the mixture weights.
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

