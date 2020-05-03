#' @title Fit Empirical Bayes Multivariate Normal Means Model
#'
#' @description This function implements an empirical Bayes method for
#' fitting a multivariate normal means model. This method is closely
#' related to approaches for multivariate density deconvolution
#' (Sarkar \emph{et al}, 2018), so it can also be viewed as a method
#' for multivariate density deconvolution, In the model, each
#' m-dimensional observation \eqn{x} is drawn from a mixture of
#' multivariate normals, \eqn{x ~ w_1 N(0, T_1) + ... + w_k N(0,
#' T_k)}, where \eqn{k} is the number of mixture components, the
#' \eqn{w_j}'s are the ixture weights, and each \eqn{T_j = S + U_j} is
#' a covariance matrix. This is the marginal density derived from a
#' model in which each \eqn{x} is multivariate normal with mean
#' \eqn{y} and covariance \eqn{S}, and the underlying, or "latent",
#' signal \eqn{y} is in turn modeled by a mixture prior in which each
#' mixture component \eqn{j} is multivariate normal with zero mean and
#' covariance matrix \eqn{S_j}. The "Extreme Deconvolution" model
#' (Bovy \emph{et al}, 2011) is a slight generalization of this
#' multivariate normal means model; it allows for the mixture
#' components to have nonzero means, and it allows one to specify a
#' linear projection of the underlying signal onto the observed
#' signal. Therefore, this method also implements a useful special
#' case of Extreme Deconvolution.
#'
#' @details The multivariate normal means model is fit by
#' expectation-maximization (EM). Two EM variants are implemented: the
#' EM algorithm described by Bovy \emph{et al} (2011); and
#' \code{"teem"} (truncated eigenvalue expectation-maximization), in
#' which the M-step update for each covariance matrix \code{U[[j]]} is
#' solved by truncating the eigenvalues in a spectral decomposition of
#' the (unconstrained) maximimum likelihood estimate (MLE) of
#' \code{U[j]]}. The latter method provides greater freedom in the
#' updates for U, and is expected to yield better fits. The choice of
#' M-step update is specified by the \code{control$update.U}
#' optimization setting; continuing reading for details.
#'
#' Using this function requires some care; currently only minimal
#' argument checking is performed. See the documentation and examples
#' for guidance.
#' 
#' @param X The n x m data matrix, in which each row is an
#'   m-dimensional observation.
#'
#' @param k An integer 2 or greater specifying the number of
#'   components in the mixture-of-normals prior. This only needs to be
#'   provided if neither \code{U} nor \code{w} are provided.
#'
#' @param w A numeric vector of length k giving initial estimates of
#'   the prior mixture weights. All entries must be non-negative, but
#'   need not sum to 1; the mixture weights are automatically normalized
#'   to sum to 1. If not provided, the mixture weights are set to
#'   uniform.
#'
#' @param U Describe input argument "U" here.
#' 
#' @param S An m x m matrix giving the initial estimate of the
#'   residual covariance matrix.
#'
#' @param control A list of parameters controlling the behaviour of
#'   the model fitting algorithm. See \sQuote{Details}.
#'
#' @return A list object with the following elements:
#'
#' \item{item1}{Describe "item1" here.}
#' 
#' \item{item2}{Describe "item2" here.}
#' 
#' @references
#'
#' J. Bovy, D. W. Hogg and S. T. Roweis (2011). Extreme Deconvolution:
#' inferring complete distribution functions from noisy, heterogeneous
#' and incomplete observations. \emph{Annals of Applied Statistics},
#' \bold{5}, 1657–1677. doi:10.1214/10-AOAS439
#'
#' A. Sarkar, D. Pati, A. Chakraborty, B. K. Mallick and R. J. Carroll
#' (2018). Bayesian semiparametric multivariate density
#' deconvolution. \emph{Journal of the American Statistical
#' Association} \bold{113}, 401–416. doi:10.1080/01621459.2016.1260467
#'
#' J. Won, J. Lim, S. Kim and B. Rajaratnam
#' (2013). Condition-number-regularized covariance estimation.
#' \emph{Journal of the Royal Statistical Society, Series B} \bold{75},
#' 427–450. doi:10.1111/j.1467-9868.2012.01049.x
#' 
#' @useDynLib mvebnm
#'
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
# Input parameters:
# U: initial estimates for covariance matrices, list of m by m matrices.
# S: T = U + S, where T is the covariance matrix for zscore_hat. And in TEEM, S = I(identity matrix)
# eps: we usually specify a small number for eps, say 1e-8, used to resolve numerical issues.
# maxiter: maximum number of iterations
# tol: criteria for convergence
#
mvebnm <- function (X, k, w, U, S = diag(ncol(X)), control = list(),
                    verbose = TRUE) {
    
  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check the input data matrix, X.
  if (!(is.numeric(X) & is.matrix(X)))
    stop("Input argument \"X\" should be a numeric matrix")

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)
  
  # Check and process input argument "k" specifying the number of
  # components in the mixture-of-normals prior.
  if (missing(k)) {
    if (missing(w) & missing(U))
      stop("At least one of \"k\", \"w\" and \"U\" should be provided")
    else if (missing(w))
      k <- length(U)
    else
      k <- length(w)
  }
  if (k < 2)
    stop("The number of prior mixture components (k) should be 2 or greater")
        
  # Check and process input argument "w" giving the initial estimates
  # of the mixture weights. Make sure the mixture weights are all
  # non-negative and sum to 1.
  if (missing(w))
    w <- rep(1,k)
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\" should be a vector of non-negative weights ",
         "of length \"k\"")
  w <- w/sum(w)
  
  # Check and process input argument "U".
  U.is.valid <- FALSE
  if (is.missing(U)) {
    # TO DO.
  }
  if (is.list(U))
    if (all(sapply(U,is.matrix))) {
      #
      # TO DO: It is important to check that the U's yield valid
      # marginal covariance matrices (the T's).
      #
      U.is.valid <- TRUE
    }
  if (!U.is.valid)
    stop("Input argument \"U\" should be list in which ...")
  
  # Check input argument "S" giving the initial estimate of the
  # residual covariance matrix.
  if (!(is.matrix(S) & nrow(S) == m & ncol(S) == m))
    stop("Input argument \"S\" should be an m x m matrix")

  # Check and process the optimization settings.
  control <- modifyList(mvebnm_control_default(),control,keep.null = TRUE)
  
  # Output the mixture weights (w), the prior covariance matrices (U),
  # and the residual covariance (S).
  for (i in 1:k) {
    rownames(U[[i]]) <- colnames(X)
    colnames(U[[i]]) <- colnames(X)
  }
  names(w)    <- names(U)
  rownames(S) <- colnames(X)
  colnames(S) <- colnames(X)
  fit         <- list(w = w,U = U,S = S)
  class(fit)  <- c("mvebnm_fit","list")
  return(fit)
  
  # initialize progress to store progress at each iteration
  progress = data.frame(iter = 1:maxiter,
                        obj  = rep(0,maxiter),
                        maxd = rep(0,maxiter))
                        
  T = c()
  
  # get a list of T. (T = U + S)
  for (i in 1:k)
    T[[i]] = U[[i]] + S

  if (verbose)
    cat("iter         objective max.diff\n")

  for (iter in 1:maxiter){
      
    # store parameters and likelihood in the previous step
    w0  <- w
    T0  <- T

    # E-step
    # ------
    # Calculate posterior probabilities using the current mu and sigmas.
    logP = matrix(0, nrow = n, ncol = k)
    for (j in 1:k){
      # log-density
      logP[,j] = log(w[j]) + dmvnorm(X, sigma = T[[j]], log = TRUE)
    }
    P = t(apply(logP, 1, function(x) normalizelogweights(x)))

    # M-step
    # ------
    # update covariance matrix with constraints
    for (j in 1:k) {
      T[[j]] = t(X)%*%(P[,j]*(X))/sum(P[,j])
      T[[j]] = shrink.cov(T[[j]], eps)
    }

    # update mixture weight
    w = colSums(P)/n

    # Compute log-likelihood at the current estimates
    f = loglik.compute(X, w, T)
    
    d = max(abs(w - w0))
    progress[iter,"obj"]  <- f
    progress[iter,"maxd"] <- d


    if (verbose)
      cat(sprintf("%4d %+0.10e %0.2e\n",iter,f,d))

    if (d < tol)
      break
  }
  
  for (i in 1:k){
    U[[i]] = T[[i]] - S
  }
  return(list(w = w, U = U, progress = progress[1:iter, ]))
}

#' @rdname mvebnm
#'
#' @export
#' 
mvebnm_control_default <- function()
  list(update.U = "teem",  # One of "em", "teem" or "none".
       update.w = "em",    # Either "em" or "none".
       update.S = "em",    # Either "em" or "none".
       maxiter  = 1000,
       eps      = 1e-15,
       tol      = 1e-6)
