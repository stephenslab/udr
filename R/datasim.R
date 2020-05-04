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
#' @importFrom mvtnorm rmvnorm
#' 
#' @export
#' 
simulate_ebmvnm_data <- function (n, w, U, S = diag()) {

  # Get the dimension of the data points (m) and the number of
  # mixture components (k).
  m <- nrow(S)
  k <- length(w)

  # Check the residual covariance matrix, S.
  if (!(is.matrix(S) & nrow(S) == m & ncol(S) == m))
    stop("Input argument \"S\" should be an m x m matrix")
  
  # Check the prior covariance matrices, U.
  U.is.valid <- FALSE
  if (is.list(U))
    if (all(sapply(U,is.matrix))) 
      U.is.valid <- verify.marginal.covariances(U,S)
  if (!U.is.valid)
    stop("Input argument \"U\" should be list in which each list element ",
         "U[[i]] is a matrix such that S + U[[i]] is symmetric positive ",
         "definite")

  # Check the mixture weights, w.
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\" should be a vector of non-negative weights ",
         "of length \"k\"")
  w <- w/sum(w)
  
  # Initialize storage for the data.
  X <- matrix(0,n,m)

  # Repeat for each sample.
  for (i in 1:n) {

    # Draw a mixture component according to the mixture weights.
    j <- sample(k,1,prob = w)
    
    # Draw the data point.
    X[i,] <- rmvnorm(1,sigma = S + U[[j]])
  }
  
  return(X)
}
