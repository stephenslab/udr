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
#' @return An n x m matrix in which each row is a draw from the
#'   multivariate normal means model.
#' 
#' @seealso \code{\link{mvebnm}}
#'
#' @importFrom mvtnorm rmvnorm
#' 
#' @export
#' 
simulate_ebmvnm_data <- function (n, w, U, S) {

  # Check the residual covariance matrix, S.
  if (!(is.matrix(S) && is.semidef(S)))
    stop("Input argument \"S\" should be a positive semi-definite matrix")
  
  # Check the prior covariance matrices, U.
  U.is.valid <- FALSE
  if (is.list(U))
    if (all(sapply(U,is.matrix))) 
      U.is.valid <- verify.prior.covariances(U,S)
  if (!U.is.valid)
    stop("Input argument \"U\" should be list in which each list element ",
         "U[[i]] is a (symmetric) positive semi-definite matrix, and",
         "S + U[[i]] is symmetric positive definite")

  # Get the dimension of the data points (m) and the number of
  # mixture components (k).
  m <- nrow(S)
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
    X[i,] <- rmvnorm(length(i),sigma = S + U[[j]])
  }
  colnames(X) <- rownames(S)
  return(X)
}
