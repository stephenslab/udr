#' @title Initialize Ultimate Deconvolution Model
#'
#' @description Initialize an Ultimate Deconvolution model. See
#'   \code{\link{ud_fit}} for background and model definition.
#' 
#' @param X The n x m data matrix, in which each row of the matrix is
#'   an m-dimensional data point. m should be 2 or greater.
#'
#' @param V Either an m x m matrix giving the initial estimate of the
#'   residual covariance matrix, or a list of m x m "standard error"
#'   matrices, one for each data point.
#'
#' @param n_rank1 A non-negative integer specifying the number of
#'   rank-1 covariance matrices included in the
#'   mixture-of-multivariate-normals prior. Initial estimates of the m x
#'   m rank-1 covariance matrices are generated at random. At most one
#'   of \code{n_rank1} and \code{U_rank1} should be provided. If neither
#'   are specified, 4 rank-1 matrices will be used.
#' 
#' @param n_unconstrained A non-negative integer specifying the number
#'   of unconstrained covariance matrices included in the
#'   mixture-of-multivariate-normals prior. Initial estimates of the m x
#'   m covariance matrices are generated at random. At most one of
#'   \code{n_unconstrained} and \code{U_unconstrained} should be
#'   provided. If neither are specified, 4 unconstrained matrices are
#'   used.
#'
#' @param U_scaled A list specifying initial estimates of the scaled
#'   covariance matrices in the mixture-of-multivariate-normals prior.
#' 
#' @param U_rank1 A list specifying initial estimates of the rank-1
#'   matrices in the mixture-of-multivariate-normals prior. At most one
#'   of \code{n_rank1} and \code{U_rank1} should be provided.
#'
#' @param U_unconstrained A list specifying initial estimates of the
#'   unconstrained matrices in the mixture-of-multivariate-normals
#'   prior. At most one of \code{n_unconstrained} and
#'   \code{U_unconstrained} should be provided.
#' 
#' @param w A numeric vector of length k giving initial estimates of
#'   the prior mixture weights, where k is equal to the number of
#'   matrices in the mixture-of-multivariate-normals prior. All entries
#'   must be non-negative, but need not sum to 1; the mixture weights
#'   are automatically normalized to sum to 1. If not provided, the
#'   mixture weights are set to uniform.
#'
#' @param version R and C++ implementations of some numerical
#'   computations are provided; these are selected with \code{version =
#'   "R"} and \code{version = "Rcpp"}.
#' 
#' @return An object of class "ud_fit". See \code{\link{ud_fit}} for
#'   details.
#'
#' @export
#'
ud_init <- function (X, V = diag(ncol(X)), n_rank1, n_unconstrained,
                     U_scaled = list(indep = diag(ncol(X)),
                                     iden = matrix(1,ncol(X),ncol(X))),
                     U_rank1, U_unconstrained, w, version = "Rcpp") {

  # Check the input data matrix, X.
  if (!(is.matrix(X) & is.numeric(X)))
    stop("Input argument \"X\" should be a numeric matrix")
  if (ncol(X) < 2)
    stop("Input argument \"X\" should have at least 2 columns")
  
  # Get the number of rows (n) and columns (m) of the data matrix.
  n <- nrow(X)
  m <- ncol(X)

  # Process the "n" and "U" input arguments (that is, n_rank1,
  # n_unconstrained, U_rank1 and U_unconstrained). First, verify that
  # at most one of n_rank1 and U_rank1 is provided, and that at most
  # one of n_unconstrained and U_unconstrained is provided.
  if (!missing(n_rank1) & !missing(U_rank1))
    stop("At most one of n_rank1 and U_rank1 should be provided")
  if (!missing(n_unconstrained) & !missing(U_unconstrained))
    stop("At most one of n_unconstrained and U_unconstrained should be ",
         "provided")
  if (missing(U_rank1)) {
    if (missing(n_rank1))
      n_rank1 <- 4
    if (n_rank1 == 0) 
      U_rank1 <- NULL
    else {
        
      # Randomly initialize the rank-1 covariance matrices.
      U_rank1 <- vector("list",n_rank1)
      for (i in 1:n_rank1)
        U_rank1[[i]] <- sim_rank1(m)
    }
  }
  if (missing(U_unconstrained)) {
    if (missing(n_unconstrained))
      n_unconstrained <- 4
    if (n_unconstrained == 0)
      U_unconstrained <- NULL
    else {
    
      # Randomly initialize the unconstrained covariance matrices.
      U_unconstrained <- vector("list",n_unconstrained)
      for (i in 1:n_unconstrained)
        U_unconstrained[[i]] <- sim_unconstrained(m)
    }
  }

  # Force all the rank-1 covariance matrices to be rank 1.
  n_scaled        <- length(U_scaled)
  n_rank1         <- length(U_rank1)
  n_unconstrained <- length(U_unconstrained)
  if (n_rank1 > 0)
    for (i in 1:n_rank1)
      U_rank1[[i]] <- getrank1(U_rank1[[i]])

  # Fill out the attributes of the covariance matrices and, if
  # necessary, label the individual matrices.
  if (n_scaled > 0) {
    if (is.null(names(U_scaled)))
      names(U_scaled) <- paste("scaled",1:n_scaled,sep = "_")
    for (i in 1:n_scaled)
      attr(U_scaled[[i]],"covtype") <- "scaled"
  }
  if (n_rank1 > 0) {
    if (is.null(names(U_rank1)))
      names(U_rank1) <- paste("rank1",1:n_rank1,sep = "_")
    for (i in 1:n_rank1)
      attr(U_rank1[[i]],"covtype") <- "rank1"
  }
  if (n_unconstrained > 0) {
    if (is.null(names(U_unconstrained)))
      names(U_unconstrained) <- paste0("unconstrained",1:n_unconstrained)
    for (i in 1:n_unconstrained)
      attr(U_unconstrained[[i]],"covtype") <- "unconstrained"
  }

  # Verify that all scaled and unconstrained matrices are
  # positive semi-definite.
  if (n_scaled > 0)
    for (i in 1:n_scaled)
      if (!issemidef(U_scaled[[i]]))
        stop("All U_scaled matrices should be positive semi-definite")
  if (n_unconstrained > 0)
    for (i in 1:n_unconstrained)
      if (!issemidef(U_unconstrained[[i]]))
        stop("All U_unconstrained matrices should be positive semi-definite")

  # Combine the prior covariances matrices into a single list.
  U <- c(U_scaled,U_rank1,U_unconstrained)
  k <- length(U)

  # Check input argument V.
  #
  # TO DO: Handle case when each data point has a separate V.
  #
  if (!issemidef(V))
    stop("V should be a positive semi-definite matrix")
  
  # Check and process input w.
  if (missing(w))
    w <- rep(1,k)
  if (length(w) != k)
    stop("Input argument \"w\" should be a numeric vector with one element ",
         "per prior covariance matrix")
   # Check the mixture weights, w.
  if (!(is.numeric(w) & all(w >= 0)))
    stop("Input argument \"w\" should be a vector containing non-negative ",
         "weights")
  w <- w/sum(w)
  names(w) <- names(U)
  
  # Add row and column names to the matrices.
  if (is.matrix(V)) {
    rownames(V) <- colnames(X)
    colnames(V) <- colnames(X)
  } else {
    for (i in 1:n) {
      rownames(V[[i]]) <- colnames(X)
      colnames(V[[i]]) <- colnames(X)
    }
  }
  for (i in 1:k) {
    rownames(U[[i]]) <- colnames(X)
    colnames(U[[i]]) <- colnames(X)
  }

  # Compute the log-likelihood.
  loglik <- loglik_ud(X,w,array(simplify2array(U),c(m,m,k)),V,version)

  # Initialize the data frame for keeping track of the algorithm's
  # progress over time.
  progress <- as.data.frame(matrix(0,0,6))
  names(progress) <- c("iter","loglik","delta.w","delta.v","delta.u","timing")
  
  # Finalize the output.
  fit <- list(V = V,U = U,w = w,loglik = loglik,progress = progress)
  class(fit) <- c("ud_fit","list")
  return(fit)
}
