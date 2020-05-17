#' @title Fit Empirical Bayes Multivariate Normal Means Model
#'
#' @description This function implements an empirical Bayes method for
#' fitting a multivariate normal means model. This method is closely
#' related to approaches for multivariate density deconvolution
#' (Sarkar \emph{et al}, 2018), so it can also be viewed as a method
#' for multivariate density deconvolution. 
#'
#' @details In the multivariate normal means model, each m-dimensional
#' observation \eqn{x} is drawn from a mixture of multivariate
#' normals, \eqn{x ~ w_1 N(0, T_1) + ... + w_k N(0, T_k)}, where
#' \eqn{k} is the number of mixture components, the \eqn{w_j}'s are
#' the mixture weights, and each \eqn{T_j = S + U_j} is a covariance
#' matrix. This is the marginal density derived from a model in which
#' each \eqn{x} is multivariate normal with mean \eqn{y} and
#' covariance \eqn{S}, and the underlying, or "latent", signal \eqn{y}
#' is in turn modeled by a mixture prior in which each mixture
#' component \eqn{j} is multivariate normal with zero mean and
#' covariance matrix \eqn{S_j}. The "Extreme Deconvolution" (ED) model
#' (Bovy \emph{et al}, 2011) is a slight generalization of this
#' multivariate normal means model; the ED model allows for the
#' mixture components to have nonzero means, sample-specific residual
#' covariances, and it allows one to specify a linear projection of
#' the underlying signal onto the observed signal. Therefore, this
#' method also implements a useful special case of Extreme
#' Deconvolution.
#'
#' The multivariate normal means model is fit by
#' expectation-maximization (EM). The \code{control} argument is a
#' list in which any of the following named components will override
#' the default optimization algorithm settings (as they are defined by
#' \code{mvebnm_control_default}):
#' 
#' \describe{
#'
#' \item{\code{update.w}}{When \code{update.w = "em"},
#' maximum-likelihood estimates of the mixture weights are computed
#' via an EM algorithm; when \code{update.w = "mixsqp"}, the mix-SQP
#' solver is used to update the mixture weights by maximizing the
#' likelihood with the other parameters fixed; when \code{update.w =
#' "none"}, the mixture weights are not updated.}
#'
#' \item{\code{update.U}}{This setting determines the algorithm used
#' to estimate the prior covariance matrices. Two EM variants are
#' implemented: \code{update.S = "em"}, the EM algorithm described by
#' Bovy \emph{et al} (2011); and \code{update.S = "teem"} (truncated
#' eigenvalue expectation-maximization), in which the M-step update
#' for each covariance matrix \code{U[[j]]} is solved by truncating
#' the eigenvalues in a spectral decomposition of the (unconstrained)
#' maximimum likelihood estimate (MLE) of \code{U[j]]}. The latter
#' method provides greater freedom in the updates for U, and is
#' expected to yield better fits, and converge more quickly to a good
#' fit.}
#'
#' \item{\code{update.S}}{When \code{update.S = "em"},
#' maximum-likelihood estimate of the residual covariance matrix is
#' computed via EM; when \code{update.S = "none", the residual
#' covariance parameter is not updated.}}
#' 
#' \item{\code{version}}{R and C++ implementations of the model
#' fitting algorithm are provided; these are selected with
#' \code{version = "R"} and \code{version = "Rcpp"}.}
#' 
#' \item{\code{maxiter}}{The upper limit on the number of EM updates
#' to perform.}
#'
#' \item{\code{tol}}{Convergence tolerance for the EM algorithm; the
#' updates are halted when the largest change in the model parameters
#' between two successive iterations of EM is less than \code{tol}.}
#'
#' \item{\code{minval}}{A small, non-negative number specifying the
#' lower bound on the eigenvalues of the prior covariance matrices
#' \code{U}.}}
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
#' @param fit0 A previous mvebnm fit. This is useful for "re-fitting"
#'   the model. When a value for this argument is given, the other
#'   inputs \code{k}, \code{w}, \code{U} and \code{S} are not needed.
#'
#' @param w A numeric vector of length k giving initial estimates of
#'   the prior mixture weights. All entries must be non-negative, but
#'   need not sum to 1; the mixture weights are automatically normalized
#'   to sum to 1. If not provided, the mixture weights are set to
#'   uniform.
#'
#' @param U A list of length k giving initial estimates of the
#'   covariance matrices in the mixture-of-multivariate-normals prior;
#'   list element \code{U[[i]]} is the m x m covariance matrix for the
#'   ith mixture component. If not provided, the initial estimates are
#'   randomly generated by computing sample covariances on random
#'   subsets of the data \code{X}.
#' 
#' @param S An m x m matrix giving the initial estimate of the
#'   residual covariance matrix.
#'
#' @param control A list of parameters controlling the behaviour of
#'   the model fitting algorithm. See \sQuote{Details}.
#'
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at each
#'   iteration. For interpretation of the columns, see the description
#'   of the \code{progress} return value.
#'
#' @return A list object with the following elements:
#'
#' \item{w}{A vector containing the estimated prior mixture
#'   weights. When \code{control$update.w = FALSE}, this will be the
#'   same as the \code{w} provided at input.}
#'
#' \item{U}{A list containing the estimated prior covariance matrices. When
#'   \code{control$update.U = FALSE}, this will be the same as the \code{U}
#'   provided at input.}
#' 
#' \item{S}{The estimated residual covariance matrix. When
#'   \code{control$update.S = FALSE}, this will be the same as the \code{S}
#'   provided at input.}
#'
#' \item{loglik}{The log-likelihood at the current settings of the
#'   model parameters, \code{w}, \code{U} and \code{S}.}
#'   
#' \item{progress}{A data frame containing detailed information about
#'   the algorithm's progress. The columns of the data frame are:
#'   "iter", the iteration number; "loglik", the log-likelihood at the
#'   current estimates of the model parameters; "delta.w", the largest
#'   change in the mixture weights; "delta.u", the largest change in the
#'   prior covariance matrices; "delta.s", the largest change in the
#'   residual covariance matrix; and "timing", the elapsed time in
#'   seconds (recorded using \code{\link{proc.time}}).}
#'
#' @examples
#' library(mvtnorm)
#' set.seed(1)
#' X   <- rmvt(1000,diag(2),df = 4)
#' fit <- mvebnm(X,k = 10,S = diag(2))
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
mvebnm <- function (X, k, fit0, w, U, S = diag(ncol(X)), control = list(),
                    verbose = TRUE) {
    
  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check the input data matrix, X.
  X <- as.matrix(X)  
  if (!is.numeric(X))
    stop("Input argument \"X\" should be a numeric matrix")

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Set up the data structure for keeping track of the algorithm's progress.
  progress <- as.data.frame(matrix(as.numeric(NA),control$maxiter,6))
  names(progress) <- c("iter","loglik","delta.w","delta.u","delta.s","timing")
  progress[,"iter"] <- 1:control$maxiter
  
  # Verify input argument "fit0", if provided, and extract the model
  # parameters and other information.
  if (!missing(fit0)) {
    if (!(missing(k) & missing(w) & missing(U) & missing(S)))
      stop("If \"fit0\" is given, do not also provide \"k\", \"w\", \"U\" ",
           "and \"S\"")
    if (!(is.list(fit0) & inherits(fit0,"mvebnm_fit")))
       stop("Input argument \"fit0\" should be an object of class ",
            "\"mvbnm_fit\", such as an output of \"mvebnm\"")
    w             <- fit0$w
    U             <- fit0$U
    S             <- fit0$S
    iter0         <- nrow(fit0$progress)
    progress$iter <- progress$iter + iter0
    progress      <- rbind(fit0$progress,progress)
  } else
    iter0 <- 0

  # Check and process input argument "k" specifying the number of
  # components in the mixture-of-normals prior.
  if (missing(k)) {
    if (missing(w) & missing(U) & missing(fit0))
      stop("At least one of \"k\", \"w\", \"U\" and \"fit0\" should be ",
           "provided")
    else if (!missing(w))
      k <- length(w)
    else
      k <- length(U)
  } else if (!(missing(w) & missing(U) & missing(fit0)))
    stop("Do not specify \"k\" if \"w\", \"U\" or \"fit0\" are already given")
  if (k < 2)
    stop("The number of prior mixture components (k) should be 2 or greater")
      
  # Check input argument "S" giving the initial estimate of the
  # residual covariance matrix.
  S <- as.matrix(S)
  if (!(nrow(S) == m & ncol(S) == m && is.semidef(S)))
    stop("Input argument \"S\" (or fit$S) should be an m x m positive ",
         "semi-definite matrix, where m is the number of columns in \"X\"")

  # Check and process input argument "U" giving the initial estimates
  # of the prior covariance matrices. If this input is not provided,
  # the initial estimates are randomly generated by computing sample
  # covariances on random subsets of the data.
  if (missing(fit0))
    if (missing(U))
      U <- generate.random.covariances(X,k)
  U <- lapply(U,as.matrix)
  if (!(all(sapply(U,is.matrix)) && verify.prior.covariances(U,S)))
    stop("Input argument \"U\" (or fit0$U) should be list in which each ",
         "list element U[[i]] is a (symmetric) positive semi-definite ",
         "matrix, and S + U[[i]] is symmetric positive definite")
  mixture.labels <- names(U)
  U <- array(simplify2array(U),c(m,m,k))  

  # Check and process input argument "w" giving the initial estimates
  # of the mixture weights. Make sure the mixture weights are all
  # non-negative and sum to 1.
  if (missing(fit0))
    if (missing(w))
       w <- rep(1,k)
  if (!(is.numeric(w) & length(w) == k & all(w >= 0)))
    stop("Input argument \"w\"(or fit0$w)  should be a vector of ",
         "length \"k\" containing non-negative weights")
  w <- w/sum(w)
  
  # Check and process the optimization settings.
  control <- modifyList(mvebnm_control_default(),control,keep.null = TRUE)

  # Give an overview of the optimization settings.
  if (verbose) {
    cat(sprintf("Fitting %d-component mvebnm to %d x %d data matrix ",k,n,m))
    cat("with these settings:\n")
    cat(sprintf("max %d updates, conv tol %0.1e ",control$maxiter,control$tol))
    cat(sprintf("(mvebnm 0.1-66, \"%s\").\n",control$version))
    cat(sprintf("updates: w (mix weights) = %s; ",control$update.w))
    cat(sprintf("U (prior cov) = %s; ",control$update.U))
    cat(sprintf("S (resid cov) = %s\n",control$update.S))
  }
  
  # RUN UPDATES
  # -----------
  if (verbose)
    cat("iter          log-likelihood |w - w'| |U - U'| |S - S'|\n")
  fit <- mvebnm_main_loop(X,w,U,S,iter0,progress,control,verbose)

  # Output the parameters of the updated model (w, U, S), the
  # log-likelihood of the updated model (loglik), and a record of the
  # algorithm's progress over time (progress).
  fit$loglik <- loglik_mvebnm(X,fit$w,fit$U,fit$S,control$version)
  fit$U      <- array2list(fit$U)
  fit$S      <- drop(fit$S)
  for (i in 1:k) {
    rownames(fit$U[[i]]) <- colnames(X)
    colnames(fit$U[[i]]) <- colnames(X)
  }
  names(fit$w)    <- mixture.labels
  names(fit$U)    <- mixture.labels
  rownames(fit$S) <- colnames(X)
  colnames(fit$S) <- colnames(X)
  class(fit)      <- c("mvebnm_fit","list")
  return(fit)
}

# This implements the core part of mvebnm.
mvebnm_main_loop <- function (X, w, U, S, iter0, progress, control, verbose) {

  # Get the number of components in the mixture prior.
  k <- length(w)

  # Iterate the EM updates.
  iterations <- iter0 + 1:control$maxiter
  for (iter in iterations) {
    t1 <- proc.time()

    # E-step
    # ------
    # Compute the n x k matrix of posterior mixture assignment
    # probabilities given current estimates of the model parameters.
    P <- compute_posterior_probs(X,w,U,S,control$version)

    # M-step
    # ------
    # Compute the M-step update for the residual covariance (S), if
    # relevant.
    if (control$update.S == "em")
      Snew <- update_resid_covariance(X,U,S,P,control$version)
    else if (control$update.S == "none")
      Snew <- S
    
    # Update the prior covariance matrices (U), if requested.
    if (control$update.U  == "em")
      Unew <- update_prior_covariance_ed(X,U,S,P,control$version)
    else if (control$update.U == "teem")
      Unew <- update_prior_covariance_teem(X,S,P,control$minval,
                                           control$version)
    else if (control$update.U == "none")
      Unew <- U
    
    # Update the mixture weights (w), if requested. Since the "mixsqp"
    # update does not use the posterior probabilities computed in the
    # E-step, it can make use of the new estimates of the other model
    # parameters.
    if (control$update.w == "em")
      wnew <- update_mixture_weights_em(P)
    else if (control$update.w == "mixsqp")
      wnew <- update_mixture_weights_mixsqp(X,Snew,Unew)
    else if (control$update.w == "none")
      wnew <- w
    
    # Update the "progress" data frame with the log-likelihood and
    # other quantities, and report the algorithm's progress to the
    # console if requested.
    loglik <- loglik_mvebnm(X,wnew,Unew,Snew,control$version)
    dw     <- max(abs(wnew - w))
    dS     <- max(abs(Snew - S))
    dU     <- max(abs(Unew - U))
    t2     <- proc.time()
    progress[iter,"loglik"]  <- loglik
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    progress[iter,"delta.s"] <- dS 
    progress[iter,"timing"]  <- t2["elapsed"] - t1["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.16e %0.2e %0.2e %0.2e\n",iter,loglik,dw,dU,dS))

    # Apply the parameter updates, and check convergencce.
    w <- wnew
    U <- Unew
    S <- Snew
    if (max(dw,dU,dS) < control$tol)
      break
  }

  # Output the parameters of the updated model (w, U, S) and a record
  # of the algorithm's progress over time ("progress").
  return(list(w = w,U = U,S = S,progress = progress[1:iter,]))
}

# Compute the n x k matrix of posterior mixture assignment
# probabilities given current estimates of the model parameters. This
# implements the "E step" in the EM algorithm for fitting the mvebnm
# model. These posterior computations are implemented in R (version =
# "R") and C++ (version = "Rcpp").
compute_posterior_probs <- function (X, w, U, S, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R")
    P <- compute_posterior_probs_helper(X,w,U,S)
  else if (version == "Rcpp")
    P <- compute_posterior_probs_rcpp(X,w,U,S)
  return(P)
}

# This implements the calculations for compute_posterior_probs.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_helper <- function (X, w, U, S) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  k <- length(w)

  # Compute the log-probabilities, stored in an n x k matrix.
  P <- matrix(0,n,k)
  for (i in 1:k)
    P[,i] = log(w[i]) + dmvnorm(X,sigma = S + U[,,i],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(softmax(P))
}

# Perform an M-step update for the mixture weights in the
# mixture-of-multivariate-normals prior. Odd things can happen if the
# mixture weights are too small, so we avoid this by adding a small
# positive scalar to all the weights.
update_mixture_weights_em <- function (P) {
  n <- nrow(P)
  w <- colSums(P)/n
  w <- pmax(w,1e-4)  
  return(w/sum(w))
}

# Update the mixture weights using the "mix-SQP" algorithm. Odd things
# can happen if the mixture weights are too small, so we avoid this by
# adding a small positive scalar to all the weights.
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom mixsqp mixsqp
#' 
update_mixture_weights_mixsqp <- function (X, S, U) {
  n <- nrow(X)
  k <- dim(U)[3]
  L <- matrix(0,n,k)
  for (i in 1:k)
    L[,i] <- dmvnorm(X,sigma = S + U[,,i],log = TRUE)
  out <- mixsqp(L,log = TRUE,control = list(verbose = FALSE))
  w   <- pmax(out$x,1e-4)
  return(w/sum(w))
}

# Perform an M-step update for the prior covariance matrices using the
# update forumla derived in Bovy et al (2011). The calculations are
# implemented in both R (version = "R") and C++ (version = "Rcpp").
update_prior_covariance_ed <- function (X, U, S, P, version = c("Rcpp","R")) {
  version <- match.arg(version)
  U0 <- U
  if (version == "R") {
    k <- ncol(P)
    for (i in 1:k)
      U[,,i] <- update_prior_covariance_ed_helper(X,U[,,i],S,P[,i])
  } else if (version == "Rcpp")
    U <- update_prior_covariances_ed_rcpp(X,U,S,P)
  return(U)
}

# Perform an M-step update for the prior covariance matrices using the
# eigenvalue-truncation technique described in Won et al (2013). The
# calculations are implemented in both R (version = "R") and C++
# (version = "Rcpp").
update_prior_covariance_teem <- function (X, S, P, minval,
                                          version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R") {
    m <- ncol(X)
    k <- ncol(P)
    U <- array(0,dim = c(m,m,k))
    for (i in 1:k)
      U[,,i] <- update_prior_covariance_teem_helper(X,S,P[,i],minval)
  } else if (version == "Rcpp")
    U <- update_prior_covariances_teem_rcpp(X,S,P,minval)
  return(U)
}

# Perform an M-step update for the residual covariance matrix, S.
update_resid_covariance <- function (X, U, S, P, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R")
    S <- update_resid_covariance_helper(X,U,S,P)
  else if (version == "Rcpp")
    S <- update_resid_covariance_rcpp(X,U,S,P)
  return(S)
}

# Perform an M-step update for the residual covariance matrix, S.
update_resid_covariance_helper <- function (X, U, S, P) {
  n    <- nrow(X)
  m    <- ncol(X)
  Snew <- matrix(0,m,m)
  for (i in 1:n) {
    out  <- compute_posterior_mvtnorm_mix(X[i,],P[i,],U,S)
    Snew <- Snew + out$S1 + tcrossprod(X[i,] - out$mu1)
  }
  return(Snew/n)
}

# Perform an M-step update for one of the prior covariance matrices
# using the update formula derived in Bovy et al (2011). Here, p is a
# vector, with one entry per row of X, giving the posterior assignment
# probabilities for the mixture component being updated.
update_prior_covariance_ed_helper <- function (X, U, S, p) {
  T <- S + U
  B <- solve(T,U)
  Y <- crossprod((sqrt(p/sum(p)) * X) %*% B)
  return(U - U %*% B + Y)
}

# Perform an M-step update for one of the prior covariance matrices
# using the eigenvalue-truncation technique described in Won et al
# (2013). Input p is a vector, with one entry per row of X, giving
# the posterior assignment probabilities for the mixture components
# being updated.
update_prior_covariance_teem_helper <- function (X, S, p, minval) {

  # Transform the data so that the residual covariance is I, then
  # compute the maximum-likelhood estimate (MLE) for T = U + I.
  R <- chol(S)
  T <- crossprod((sqrt(p/sum(p))*X) %*% solve(R))
  
  # Find U maximizing the expected complete log-likelihood subject to
  # U being positive definite. This update for U is based on the fact
  # that the covariance matrix that minimizes the likelihood subject
  # to the constraint that U is positive definite is obtained by
  # truncating the eigenvalues of T = U + I less than 1 to be 1; see
  # Won et al (2013), p. 434, the sentence just after equation (16).
  U <- shrink.cov(T,minval)

  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}

# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance S, and z is drawn from a mixture of multivariate
# normals, each with zero mean, covariance V[,,i] and weight w[i].
# Return the posterior mean (mu1) and covariance (S1) of z. Note that
# input w1 must be the vector of *posterior* mixture weights (see
# compute_posterior_probs).
compute_posterior_mvtnorm_mix <- function (x, w1, V, S) {
  m   <- length(x)
  k   <- length(w1)
  mu1 <- rep(0,m)
  S1  <- matrix(0,m,m)
  for (i in 1:k) {
    out <- compute_posterior_mvtnorm(x,V[,,i],S)
    mu1 <- mu1 + w1[i] * out$mu
    S1  <- S1  + w1[i] * (out$S + tcrossprod(out$mu))
  }
  S1 <- S1 - tcrossprod(mu1)
  return(list(mu1 = mu1,S1 = S1))
}

# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance S, and z is drawn from a multivariate normal
# distribution with mean zero and covariance V. Return the posterior
# mean (mu1) and covariance (S1) of z. These calculations will only
# work if S is positive definite.
compute_posterior_mvtnorm <- function (x, V, S) {
  m   <- length(x)
  S1  <- solve(V %*% solve(S) + diag(m)) %*% V
  mu1 <- drop(S1 %*% solve(S,x))
  return(list(mu1 = mu1,S1 = S1))
}

#' @rdname mvebnm
#'
#' @export
#' 
mvebnm_control_default <- function()
  list(update.w = "em",    # One or "em", "mixsqp" or "none".
       update.U = "teem",  # One of "em", "teem" or "none".
       update.S = "none",  # One of "em" or "none".
       version  = "Rcpp",  # One of "R" or "Rcpp".
       maxiter  = 100,
       minval   = 1e-8,
       tol      = 1e-6)
