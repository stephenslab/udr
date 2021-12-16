#' @title Initialize Ultimate Deconvolution Model Fit
#'
#' @description Initialize an Ultimate Deconvolution model fit. See
#'   \code{\link{ud_fit}} for background and model definition.
#' 
#' @param X An n x m data matrix, in which each row of the matrix is
#'   an m-dimensional data point. \code{X} should have at least 2 rows
#'   and 2 columns.
#'
#' @param V Either an m x m matrix giving the
#'   residual covariance matrix (assumed equal for every obsevation), or a list of m x m covariance matrices
#'   of length n.
#'
#' @param n_rank1 A non-negative integer specifying the number of
#'   rank-1 covariance matrices included in the mixture prior. Initial
#'   estimates of the m x m rank-1 covariance matrices are generated at
#'   random. At most one of \code{n_rank1} and \code{U_rank1} should be
#'   provided. If neither are specified, 4 random rank-1 matrices will be
#'   included.
#' 
#' @param n_unconstrained A non-negative integer specifying the number
#'   of unconstrained covariance matrices included in the mixture
#'   prior. Initial estimates of the m x m covariance matrices are
#'   generated at random. At most one of \code{n_unconstrained} and
#'   \code{U_unconstrained} should be provided. If neither are
#'   specified, 4 random unconstrained matrices will be included.
#'
#' @param U_scaled A list specifying initial estimates of the scaled
#'   covariance matrices in the mixture prior. (The defaults provide two commonly-used covariance matrices here)
#' 
#' @param U_rank1 A list specifying initial estimates of the rank-1
#'   matrices in the mixture prior. At most one of \code{n_rank1} and
#'   \code{U_rank1} should be provided. If \code{U_rank1} is not given,
#'   the rank-1 covariates are initialized at random.
#'
#' @param U_unconstrained A list specifying initial estimates of the
#'   unconstrained matrices in the mixture prior. At most one of
#'   \code{n_unconstrained} and \code{U_unconstrained} should be
#'   provided. If \code{U_unconstrained} is not given, the matrices are
#'   initialized at random.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the model initialization. See \code{\link{ud_fit}} for details.
#'
#' @return An Ultimate Deconvolution model fit. See
#'   \code{\link{ud_fit}} for details.
#'
#' @seealso \code{\link{ud_fit}}
#' 
#' @export
#'
ud_init <- function (X, V = diag(ncol(X)), n_rank1, n_unconstrained,
                     U_scaled = list(indep = diag(ncol(X)),
                                     equal = matrix(1,ncol(X),ncol(X)) +
                                                    1e-8 * diag(ncol(X))),
                     U_rank1, U_unconstrained, control = list()) {

  # Check the input data matrix, "X".
  if (!(is.matrix(X) & is.numeric(X)))
    stop("Input argument \"X\" should be a numeric matrix")
  n <- nrow(X)
  m <- ncol(X)
  if (n < 2 | m < 2)
    stop("Input argument \"X\" should have at least 2 columns and ",
         "at least 2 rows")
  
  # Check and process the optimization settings.
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
  
  # Check input argument "V".
  modified <- FALSE
  if (is.matrix(V)) {
    if (!issemidef(V)) {
      V <- makeposdef(V)
      modified <- TRUE
    }
  } else {
    if (length(V) != n)
      stop("Input argument \"V\" should either be a positive ",
           "semi-definite matrix, or a list of positive semi-definite ",
           "matrices, with one matrix per row of \"X\"")
    for (i in 1:n)
      if (!issemidef(V[[i]])) {
         V[[i]] <- makeposdef(V[[i]])
         modified <- TRUE
      }
  }
  if (modified)
    warning("Input argument \"V\" was modified because one or more ",
            "matrices were not positive semi-definite")
  
  # Process inputs n_rank1 and U_rank1. If U_rank1 is not provided,
  # randomly initialize the rank-1 covariance matrices.
  if (!missing(n_rank1) & !missing(U_rank1))
    stop("At most one of n_rank1 and U_rank1 should be provided")
  if (missing(U_rank1)) {
    if (missing(n_rank1))
      n_rank1 <- 4
    if (n_rank1 == 0) 
      U_rank1 <- NULL
    else {
      U_rank1 <- vector("list",n_rank1)
      for (i in 1:n_rank1)
        U_rank1[[i]] <- sim_rank1(m)
    }
  }

  # Process inputs n_unconstrained and U_unconstrained. If
  # U_unconstrained is not provided, randomly initialize the
  # unconstrained covariance matrices.
  if (!missing(n_unconstrained) & !missing(U_unconstrained))
    stop("At most one of n_unconstrained and U_unconstrained should be ",
         "provided")
  if (missing(U_unconstrained)) {
    if (missing(n_unconstrained))
      n_unconstrained <- 4
    if (n_unconstrained == 0)
      U_unconstrained <- NULL
    else {
      U_unconstrained <- vector("list",n_unconstrained)
      for (i in 1:n_unconstrained)
        U_unconstrained[[i]] <- sim_unconstrained(m)
    }
  }

  # Get the number of prior covariance matrices of each type.
  n_scaled        <- length(U_scaled)
  n_rank1         <- length(U_rank1)
  n_unconstrained <- length(U_unconstrained)

  # Verify that all scaled and unconstrained matrices are
  # positive semi-definite.
  modified <- FALSE
  if (n_scaled > 0)
    for (i in 1:n_scaled)
      if (!issemidef(U_scaled[[i]])) {
        U_scaled[[i]] <- makeposdef(U_scaled[[i]])
        modified <- TRUE
      }
  if (modified)
    warning("Input argument \"U_scaled\" was modified because one or more ",
            "matrices were not positive semi-definite")
  
  modified <- FALSE
  if (n_unconstrained > 0)
    for (i in 1:n_unconstrained)
      if (!issemidef(U_unconstrained[[i]])) {
        U_unconstrained[[i]] <- makeposdef(U_unconstrained[[i]])
        modified <- TRUE
      }
  if (modified)
    warning("Input argument \"U_unconstrained\" was modified because one ",
            "or more matrices were not positive semi-definite")

  # Set up the data structure for the scaled covariance matrices.
  if (n_scaled > 0) {
    if (is.null(names(U_scaled)))
      names(U_scaled) <- paste("scaled",1:n_scaled,sep = "_")
    for (i in 1:n_scaled)
      U_scaled[[i]] <- create_prior_covariance_struct_scaled(X,U_scaled[[i]])
  }

  # Set up the data structure for the rank-1 covariance matrices.
  if (n_rank1 > 0) {
    if (is.null(names(U_rank1)))
      names(U_rank1) <- paste("rank1",1:n_rank1,sep = "_")
    for (i in 1:n_rank1)
      U_rank1[[i]] <- create_prior_covariance_struct_rank1(X,U_rank1[[i]])
  }

  # Set up the data structure for the unconstrained covariance matrices.
  if (n_unconstrained > 0) {
    if (is.null(names(U_unconstrained)))
      names(U_unconstrained) <- paste0("unconstrained",1:n_unconstrained)
    for (i in 1:n_unconstrained)
      U_unconstrained[[i]] <-
        create_prior_covariance_struct_unconstrained(X,U_unconstrained[[i]])
  }

  # Combine the prior covariances matrices into a single list.
  U <- c(U_scaled,U_rank1,U_unconstrained)
  k <- length(U)
  if (k < 2)
    stop("The total number of prior covariances should be at least 2")

  # Initialize the mixture weights.
  w        <- rep(1,k)/k
  names(w) <- names(U)
  
  # Add row and column names to V.
  if (is.matrix(V)) {
    rownames(V) <- colnames(X)
    colnames(V) <- colnames(X)
  } else {
    names(V) <- rownames(X)
    for (i in 1:n) {
      rownames(V[[i]]) <- colnames(X)
      colnames(V[[i]]) <- colnames(X)
    }
  }

  # Initialize the data frame for keeping track of the algorithm's
  # progress over time.
  progress        <- as.data.frame(matrix(0,0,6))
  names(progress) <- c("iter","loglik","delta.w","delta.v","delta.u","timing")
  
  # Compute the log-likelihood and the responsibilities matrix (P), and
  # finalize the output.
  fit        <- list(X = X,V = V,U = U,w = w,
                     progress = progress)
  class(fit) <- c("ud_fit","list")
  fit <- compute_loglik_matrix(fit,control$version)
  fit <- compute_posterior_probs(fit)
  fit <- compute_loglik(fit)
  return(fit)
}
