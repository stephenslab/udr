#' @title Initialize Ultimate Deconvolution Model
#'
#' @description Initialize an Ultimate Deconvolution model fit. See
#'   \code{\link{ud_fit}} for background and model definition.
#' 
#' @param dat An n x m data matrix, in which each row of the matrix is
#'   an m-dimensional data point, or a \dQuote{mash} object, for example
#'   created by \code{\link[mashr]{mash_set_data}}. When \code{X} is a
#'   matrix it should have at least 2 rows and at least 2 columns.
#'   
#' @param V Either an m x m matrix giving the residual covariance
#'   matrix for all n data points, or a list of m x m covariance
#'   matrices of length n.
#'
#' @param n_rank1 A non-negative integer specifying the number of
#'   rank-1 covariance matrices included in the mixture prior. Initial
#'   estimates of the m x m rank-1 covariance matrices are generated at
#'   random. At most one of \code{n_rank1} and \code{U_rank1} should be
#'   provided. If neither are specified, rank-1 matrices will not be
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
#'   covariance matrices in the mixture prior. (The default setting
#'   specifies two commonly used covariance matrices.) In the case of a
#'   single scaled matrix, \code{U_scaled} may be a matrix instead of a
#'   list.
#' 
#' @param U_rank1 A list specifying initial estimates of the rank-1
#'   matrices in the mixture prior. At most one of \code{n_rank1} and
#'   \code{U_rank1} should be provided. If \code{U_rank1} is not given,
#'   the rank-1 covariates are initialized at random. In the case of a
#'   single rank-1 matrix, \code{U_rank1} may be a matrix instead of
#'   list.
#'
#' @param U_unconstrained A list specifying initial estimates of the
#'   unconstrained matrices in the mixture prior. At most one of
#'   \code{n_unconstrained} and \code{U_unconstrained} should be
#'   provided. If \code{U_unconstrained} is not given, the matrices are
#'   initialized at random. In the case of a single unconstrained
#'   matrix, \code{U_unconstrained} may be a matrix instead of a list.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the model initialization. See \code{\link{ud_fit}} for details.
#'
#' @return An Ultimate Deconvolution model fit. See
#'   \code{\link{ud_fit}} for details.
#'
#' @seealso \code{\link{ud_fit}}
#'
#' @examples
#' # Here we illustrate the different ways to initialize a UD model fit,
#' # with and without mashr. In the very simplest invocation, we supply a
#' # data matrix X.
#' set.seed(1)
#' X <- matrix(rnorm(2000),200,10)
#' fit <- ud_init(X)
#'
#' # For more sophisticated usage, first create a mash data object with
#' # mash_set_data, then call ud_init. By setting alpha = 1, we invoke
#' # the EZ ("exchangeable z-scores") model, and with the model the
#' # residual covariance matrix is the same for all effects (rows of
#' # Bhat), and is the same as the V in the mash object.
#' library(mashr)
#' V <- matrix(0.1,10,10)
#' diag(V) <- 1
#' Shat <- sqrt(matrix(rexp(2000),200,10))
#' out <- simple_sims(200,10,Shat)
#' dat <- mash_set_data(out$Bhat,out$Shat,V = V,alpha = 1)
#' fit <- ud_init(dat)
#' all.equal(dat$V,fit$V,check.attributes = FALSE)
#'
#' # By setting alpha = 0, we instead use the EE ("exchangeable effects")
#' # model. Now the residual covariance is different for each effect, and
#' # these covariances are now stored in a list with one element for
#' # every row of Bhat (or X).
#' dat <- mash_set_data(out$Bhat,out$Shat,V = V,alpha = 0)
#' fit <- ud_init(dat)
#' class(fit$V)
#' 
#' # Again we use the EE model, but since the standard errors are the
#' # same for all the effects (and the s.e.'s are equal to 1), the
#' # residual covariance matrix is now the same for all effects.
#' dat <- mash_set_data(out$Bhat,Shat = 1,V = V,alpha = 0)
#' fit <- ud_init(dat)
#' all.equal(dat$V,fit$V,check.attributes = FALSE)
#'
#' # In this final example, we customize the mixture prior to include
#' # different numbers of unconstrained and rank-1 matrices.
#' fit <- ud_init(dat,n_rank1 = 2,n_unconstrained = 4)
#' names(fit$U)
#' 
#' @export
#'
ud_init <- function (dat, V, n_rank1, n_unconstrained,
                     U_scaled = list(indep = diag(ncol(X)),
                                     equal = matrix(1,ncol(X),ncol(X)) +
                                             1e-4 * diag(ncol(X))),
                     U_rank1, U_unconstrained, control = list()) {

  # Check and process input arguments "dat" and "V".
  if (inherits(dat,"mash")) {
    if (!missing(V))
      warning("Input argument V is ignored because it is already supplied ",
              "by the \"dat\" argument")
    if (dat$alpha == 1) {

      # Use the "Exchangeable Z-scores" (EZ) model: the data matrix is
      # set to the z-scores, and V is common to all rows of X.
      X <- dat$Bhat/dat$Shat
      V <- dat$V  
    } else if (dat$alpha == 0) {

      # Use the "Exchangeable Effects" (EE) model: the data matrix is
      # set to the effect estimates (Bhat), and there is (potentially)
      # a different V for each row of X.
      X <- dat$Bhat  
      n <- nrow(X)
      if (dat$commonV) {
        
        # Check if shat is the same for all rows.
        shat <- dat$Shat[1,]
        if (max(abs(shat - dat$Shat)) < 1e-15)

          # All the rows of X share the same V.
          V <- diag(shat) %*% dat$V %*% diag(shat)
        else {

          # Each row of X has a different V.
          V <- vector("list",n)
          for (i in 1:n)
            V[[i]] = diag(dat$Shat[i,]) %*% dat$V %*% diag(dat$Shat[i,])
        }
      } else {

        # We have a different V for each row of X.
        V <- vector("list",n)
        for (i in 1:n)
          V[[i]] = diag(dat$Shat[i,]) %*% dat$V[,,i] %*% diag(dat$Shat[i,])
      }
    } else
      stop("Invalid value of \"alpha\" in mash object \"dat\"; should be ",
           "0 or 1 (see mashr function mash_set_data for more information)")
  } else if (is.matrix(dat) & is.numeric(dat)) {
    X <- dat
    if (missing(V))
      V <- diag(ncol(X))
  } else
    stop("Input argument \"dat\" should be a \"mash\" object or a ",
         "numeric matrix")
      
  # Perform some more checks of X.
  if (!(is.matrix(X) & is.numeric(X)))
    stop("X should be a numeric matrix")
  if (any(is.na(X)) | any(is.infinite(X)))
    stop("X should not have missing or infinite values")
  n <- nrow(X)
  m <- ncol(X)
  if (n < 2 | m < 2)
    stop("X should have at least 2 columns and at least 2 rows")
  
  # Check and process the optimization settings.
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)

  # Perform some more checks of V.
  modified <- FALSE
  if (is.matrix(V)) {
    if (any(is.na(V)) | any(is.infinite(V)))
      stop("V should not contain missing or infinite values")
    if (!issemidef(V)) {
      V <- makeposdef(V)
      modified <- TRUE
    }
  } else {
    if (length(V) != n)
      stop("V should either be a positive semi-definite matrix, or a list ",
           "of positive semi-definite matrices, with one matrix per row of X")
    for (i in 1:n)
      if (any(is.na(V[[i]])) | any(is.infinite(V[[i]])))
        stop("V should not contain missing or infinite values")
    if (!issemidef(V[[i]])) {
      V[[i]] <- makeposdef(V[[i]])
      modified <- TRUE
    }
  }
  if (modified)
    warning("V was modified because one or more matrices were not ",
            "positive semi-definite")
  
  # Process inputs n_rank1 and U_rank1. If U_rank1 is not provided,
  # randomly initialize the rank-1 covariance matrices.
  if (!missing(n_rank1) & !missing(U_rank1))
    stop("At most one of n_rank1 and U_rank1 should be provided")
  if (missing(U_rank1)) {
    if (missing(n_rank1))
      n_rank1 <- 0
    if (n_rank1 == 0) 
      U_rank1 <- NULL
    else {
      U_rank1 <- vector("list",n_rank1)
      for (i in 1:n_rank1)
        U_rank1[[i]] <- sim_rank1(m)
    }
  }
  if (is.matrix(U_rank1))
    U_rank1 <- list(U_rank1)
  
  # Process inputs n_unconstrained and U_unconstrained. If
  # U_unconstrained is not provided, randomly initialize the
  # unconstrained covariance matrices.
  if (!missing(n_unconstrained) & !missing(U_unconstrained))
    stop("At most one of n_unconstrained and U_unconstrained should be ",
         "provided")
  if (missing(U_unconstrained)) {
    if (missing(n_unconstrained))
      n_unconstrained <- 8
    if (n_unconstrained == 0)
      U_unconstrained <- NULL
    else {
      U_unconstrained <- vector("list",n_unconstrained)
      for (i in 1:n_unconstrained)
        U_unconstrained[[i]] <- sim_unconstrained(m)
    }
  }
  if (is.matrix(U_unconstrained))
    U_unconstrained <- list(U_unconstrained)
  
  # Process input U_scaled.
  if (is.matrix(U_scaled))
    U_scaled <- list(U_scaled)
  
  # Get the number of prior covariance matrices of each type.
  n_scaled        <- length(U_scaled)
  n_rank1         <- length(U_rank1)
  n_unconstrained <- length(U_unconstrained)
  
  # Verify there are no missing or infinite values in U_rank1.
  if (n_rank1 != 0)
    for (i in 1:n_rank1)
      if (any(is.na(U_rank1[[i]])) | any(is.infinite(U_rank1[[i]])))
        stop("U_rank1 should not contain missing or infinite values")
  
  # Verify that all scaled and unconstrained matrices are
  # positive semi-definite.
  modified <- FALSE
  if (n_scaled > 0) {
    for (i in 1:n_scaled) {
      if (any(is.na(U_scaled[[i]])) | any(is.infinite(U_scaled[[i]])))
        stop("U_scaled should not contain missing or infinite values")
      if (!issemidef(U_scaled[[i]])) {
        U_scaled[[i]] <- makeposdef(U_scaled[[i]])
        modified <- TRUE
      } 
    }
    if (modified)
      warning("U_scaled was modified because one or more matrices were not ",
              "positive semi-definite")
  }
  
  modified <- FALSE
  if (n_unconstrained > 0){
    for (i in 1:n_unconstrained) {
      if (any(is.na(U_unconstrained[[i]])) |
          any(is.infinite(U_unconstrained[[i]])))
        stop("U_unconstrained should not contain missing or infinite values")
      if (!issemidef(U_unconstrained[[i]])) {
        U_unconstrained[[i]] <- makeposdef(U_unconstrained[[i]])
        modified <- TRUE
      }
    }
    if (modified)
      warning("U_unconstrained was modified because one or more matrices ",
              "were not positive semi-definite")
  }
  
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
  progress        <- as.data.frame(matrix(0,0,7))
  names(progress) <- c("iter","loglik","loglik.pen","delta.w","delta.v",
                       "delta.u","timing")
  
  # Compute the log-likelihood and the responsibilities matrix (P), and
  # finalize the output.
  loglik     <- loglik_ud(X,w,U,V,control$version)
  fit        <- list(X = X,V = V,U = U,w = w,loglik = loglik, 
                     loglik_penalized = -Inf,logplt = as.numeric(NA),
                     progress = progress)
  class(fit) <- c("ud_fit","list")
  fit        <- compute_posterior_probs(fit,control$version)
  return(fit)
}
