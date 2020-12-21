#' @rdname ud_fit
#' 
#' @title Fit Ultimate Deconvolution Model
#'
#' @description This function implements "Ultimate Deconvolution", an
#' empirical Bayes method for fitting a multivariate normal means
#' model. This method is closely related to approaches for
#' multivariate density deconvolution (Sarkar \emph{et al}, 2018), so
#' it can also be viewed as a method for multivariate density
#' deconvolution.
#'
#' @details In the Ultimate Deconvolution (UD) model, the
#' m-dimensional observation \eqn{x} is drawn from a mixture of
#' multivariate normals, \eqn{x ~ w_1 N(0, T_1) + ... + w_k N(0,
#' T_k)}, where \eqn{k} is the number of mixture components, the
#' \eqn{w_j}'s are the mixture weights, and each \eqn{T_j = V + U_j}
#' is a covariance matrix. This is the marginal density derived from a
#' model in which \eqn{x} is multivariate normal with mean \eqn{y} and
#' covariance \eqn{V}, and the underlying, or "latent", signal \eqn{y}
#' is in turn modeled by a mixture prior in which each mixture
#' component \eqn{j} is multivariate normal with zero mean and
#' covariance matrix \eqn{U_j}. This model is a useful special case of
#' the "Extreme Deconvolution" (ED) model (Bovy \emph{et al}, 2011).
#'
#' Two variants of the UD model are implemented: one in which the
#' residual covariance \code{V} is the same for all data points, and
#' another in which V is different for each data point. In the first
#' case, this covariance matrix may be estimated.
#' 
#' The UD model is fit by expectation-maximization (EM). The
#' \code{control} argument is a list in which any of the following
#' named components will override the default optimization algorithm
#' settings (as they are defined by \code{ud_fit_control_default}):
#' 
#' \describe{
#'
#' \item{\code{weights.update}}{When \code{weights.update = "em"}, the
#' mixture weights are updated via EM; when \code{weights.update =
#' "none"}, the mixture weights are not updated.}
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
#' between two successive updates is less than \code{tol}.}
#'
#' \item{\code{update.U}}{This setting determines the algorithm used
#' to estimate the prior covariance matrices. Two EM variants are
#' implemented: \code{update.U = "ed"}, the EM algorithm described by
#' Bovy \emph{et al} (2011); and \code{update.U = "teem"} (truncated
#' eigenvalue expectation-maximization), in which the M-step update
#' for each covariance matrix \code{U[[j]]} is solved by truncating
#' the eigenvalues in a spectral decomposition of the unconstrained
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
#' \item{\code{minval}}{Minimum eigenvalue allowed in the residual
#'   covariance(s) \code{V} and the prior covariance matrices
#'   \code{U}.}}
#'
#' Using this function requires some care; currently only minimal
#' argument checking is performed. See the documentation and examples
#' for guidance.
#'
#' @param fit0 A previous Ultimate Deconvolution model fit. Typically,
#'   this will be an output from \code{\link{ud_init}} or from a
#'   previous call to \code{ud_fit}.
#'
#' @param X Describe input argument "X" here.
#' 
#' @param control A list of parameters controlling the behaviour of
#'   the model fitting and initialization. See \sQuote{Details}.
#'
#' @param verbose When \code{verbose = TRUE}, information about the
#'   algorithm's progress is printed to the console at each
#'   iteration. For interpretation of the columns, see the description
#'   of the \code{progress} return value.
#'
#' @return A list object with the following elements:
#'
#' \item{w}{A vector containing the estimated prior mixture
#'   weights. When \code{control$update.w = "none"}, this will be the
#'   same as the \code{w} provided at input.}
#'
#' \item{U}{A list containing the estimated prior covariance matrices. When
#'   \code{control$update.U = "none"}, this will be the same as the \code{U}
#'   provided at input.}
#' 
#' \item{V}{The estimated residual covariance matrix. When
#'   \code{control$update.S = "none"}, this will be the same as the \code{S}
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
#' # TO DO.
#' 
#' @references
#'
#' J. Bovy, D. W. Hogg and S. T. Roweis (2011). Extreme Deconvolution:
#' inferring complete distribution functions from noisy, heterogeneous
#' and incomplete observations. \emph{Annals of Applied Statistics},
#' \bold{5}, 1657–1677. doi:10.1214/10-AOAS439
#'
#' A. Sarkar, D. Pati, A. Chakraborty, B. K. Mallick and R. J. Carroll
#' (2018). Bayesian semiparametric multivariate density deconvolution.
#' \emph{Journal of the American Statistical Association} \bold{113},
#' 401–416. doi:10.1080/01621459.2016.1260467
#'
#' J. Won, J. Lim, S. Kim and B. Rajaratnam (2013).
#' Condition-number-regularized covariance estimation. \emph{Journal
#' of the Royal Statistical Society, Series B} \bold{75},
#' 427–450. doi:10.1111/j.1467-9868.2012.01049.x
#' 
#' @useDynLib udr
#'
#' @importFrom utils modifyList
#' @importFrom Rcpp evalCpp
#' 
#' @export
#' 
ud_fit <- function (fit0, X, control = list(), verbose = TRUE) {
    
  # CHECK & PROCESS INPUTS
  # ----------------------
  # Check input argument "fit0".
  if (!(is.list(fit0) & inherits(fit0,"ud_fit")))
    stop("Input argument \"fit0\" should be an object of class \"ud_fit\",",
         "such as an output of ud_init")
  fit <- fit0
  
  # Check the input data matrix, X.
  if (missing(X))
    X <- fit$X
  if (!(is.matrix(X) & is.numeric(X)))
    stop("Input argument \"X\" should be a numeric matrix")

  # Get the number of rows (n) and columns (m) of the data matrix,
  n <- nrow(X)
  m <- ncol(X)

  # Check and process the optimization settings.
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
  
  # Give an overview of the model fitting.
  if (verbose) {
    covtypes <- sapply(fit$U,function (x) attr(x,"covtype"))
    cat(sprintf("Performing Ultimate Deconvolution on %d x %d matrix ",n,m))
    cat(sprintf("(udr 0.3-11, \"%s\"):\n",control$version))
    if (is.matrix(fit$V))
      cat("data points are i.i.d. (same V)\n")
    else
      cat("data points are not i.i.d. (different Vs)\n")
    cat(sprintf("prior covariances: %d scaled, %d rank-1, %d unconstrained\n",
                sum(covtypes == "scaled"),
                sum(covtypes == "rank1"),
                sum(covtypes == "unconstrained")))
    cat(sprintf("mixture weights update: %s\n",control$weights.update))
    if (is.matrix(fit$V))
      cat(sprintf("residual covariance update: %s\n",control$resid.update))
    cat(sprintf("max %d updates, conv tol %0.1e\n",
                control$maxiter,control$tol))
    # cat(sprintf("U (prior cov) = %s; ",control$update.U))
  }
  
  # RUN UPDATES
  # -----------
  if (verbose)
    cat("iter          log-likelihood |w - w'| |U - U'| |V - V'|\n")
  if (is.list(fit$V))
    fit$V <- simplify2array(fit$V)
  fit$U <- simplify2array(fit$U)
  fit <- ud_fit_main_loop(X,fit$w,fit$U,fit$V,control,verbose)
  
  # Output the parameters of the updated model (w, U, V), the
  # log-likelihood of the updated model (loglik), and a record of the
  # algorithm's progress over time (progress).
  fit$progress <- rbind(fit0$progress,fit$progress)
  fit$loglik   <- loglik_ud(X,fit$w,fit$U,fit$V,control$version)
  fit$X        <- X
  fit$U        <- array2list(fit$U)
  if (!is.matrix(V))
    fit$V <- array2list(fit$V)
  # for (i in 1:k) {
  #   rownames(fit$U[[i]]) <- colnames(X)
  #   colnames(fit$U[[i]]) <- colnames(X)
  # }
  # names(fit$w)    <- mixture.labels
  # names(fit$U)    <- mixture.labels
  # rownames(fit$S) <- colnames(X)
  # colnames(fit$S) <- colnames(X)
  class(fit) <- c("ud_fit","list")
  return(fit)
}

# This implements the core part of ud_fit.
ud_fit_main_loop <- function (X, w, U, V, control, verbose) {

  # Get the number of components in the mixture prior.
  k <- length(w)

  # Set up data structures used in the loop below.
  progress <- as.data.frame(matrix(0,control$maxiter,6))
  names(progress) <- c("iter","loglik","delta.w","delta.v","delta.u","timing")

  # Iterate the EM updates.
  for (iter in 1:control$maxiter) {
    t1 <- proc.time()

    # E-step
    # ------
    # Compute the n x k matrix of posterior mixture assignment
    # probabilities given current estimates of the model parameters.
    P <- compute_posterior_probs(X,w,U,V,control$version)

    # M-step
    # ------
    # Update the residual covariance matrix.
    if (is.matrix(V)) {
      if (control$resid.update == "em")
        Vnew <- update_resid_covariance(X,U,V,P,control$version)
      else if (control$resid.update == "none")
        Vnew <- V
    }
    
    # Update the prior covariance matrices (U), if requested.
    # if (control$update.U == "ed")
    #   Unew <- update_prior_covariance_ed(X,U,S,P,control$version)
    # else if (control$update.U == "teem")
    #   Unew <- update_prior_covariance_teem(X,S,P,control$minval,
    #                                        control$version)
    # else if (control$update.U == "none")
    #   Unew <- U
    Unew <- U
    
    # Update the mixture weights (w), if requested. Since the "mixsqp"
    # update does not use the posterior probabilities computed in the
    # E-step, it can make use of the new estimates of the other model
    # parameters.
    if (control$weights.update == "em")
      wnew <- update_mixture_weights_em(P)
    else if (control$weights.update == "none")
      wnew <- w
    
    # Update the "progress" data frame with the log-likelihood and
    # other quantities, and report the algorithm's progress to the
    # console if requested.
    loglik <- loglik_ud(X,wnew,Unew,Vnew,control$version)
    dw <- max(abs(wnew - w))
    dU <- max(abs(Unew - U))
    dV <- max(abs(Vnew - V))
    t2 <- proc.time()
    progress[iter,"loglik"]  <- loglik
    progress[iter,"delta.w"] <- dw 
    progress[iter,"delta.u"] <- dU 
    progress[iter,"delta.v"] <- dV 
    progress[iter,"timing"]  <- t2["elapsed"] - t1["elapsed"]
    if (verbose)
      cat(sprintf("%4d %+0.16e %0.2e %0.2e %0.2e\n",iter,loglik,dw,dU,dV))

    # Apply the parameter updates, and check convergencce.
    w <- wnew
    U <- Unew
    V <- Vnew
    if (max(dw,dU,dV) < control$tol)
      break
  }

  # Output the parameters of the updated model, and a record of the
  # algorithm's progress over time.
  return(list(w = w,U = U,V = V,progress = progress[1:iter,]))
}

# Compute the n x k matrix of posterior mixture assignment
# probabilities given current estimates of the model parameters. This
# implements the E step in the EM algorithm for fitting the Ultimate
# Deconvolution model.
compute_posterior_probs <- function (X, w, U, V, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (is.matrix(V)) {

    # Perform the computations for the special case when the same
    # residual variance is used for all samples.
    if (version == "R")
      P <- compute_posterior_probs_iid_helper(X,w,U,V)
    else if (version == "Rcpp")
      P <- compute_posterior_probs_iid_rcpp(X,w,U,V)
  } else {

    # Perform the computations for the more general case when the
    # residual variance is not necessarily the same for all samples.
    if (version == "R")
      P <- compute_posterior_probs_general_helper(X,w,U,V)
    else if (version == "Rcpp")
      P <- compute_posterior_probs_general_rcpp(X,w,U,V)
  }
  return(P)
}

# This implements the calculations for compute_posterior_probs for the
# special case when the same residual covariance matrix is used for
# all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_iid_helper <- function (X, w, U, V) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  k <- length(w)

  # Compute the log-probabilities, stored in an n x k matrix.
  P <- matrix(0,n,k)
  for (j in 1:k)
    P[,j] = log(w[j]) + dmvnorm(X,sigma = V + U[,,j],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(softmax(P))
}

# This implements the calculations for compute_posterior_probs for the
# more general case when the residual covariances are not necessaily
# the same for all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_general_helper <- function (X, w, U, V) {
      
  # Get the number of samples (n) and the number of components in the
  # mixture prior (k).
  n <- nrow(X)
  k <- length(w)

  # Compute the log-probabilities, stored in an n x k matrix.
  P <- matrix(0,n,k)
  for (i in 1:n)
    for (j in 1:k)
      P[i,j] = log(w[j]) + dmvnorm(X[i,],sigma = V[,,i] + U[,,j],log = TRUE)

  # Normalize the probabilities so that each row of P sums to 1.
  return(softmax(P))
}

# Perform an M-step update for the mixture weights in the
# mixture-of-multivariate-normals prior.
update_mixture_weights_em <- function (P)
  colSums(P)/nrow(P)

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

# Perform an M-step update for the residual covariance matrix.
update_resid_covariance <- function (X, U, V, P, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "R")
    V <- update_resid_covariance_helper(X,U,V,P)
  else if (version == "Rcpp")
    V <- update_resid_covariance_rcpp(X,U,V,P)
  return(V)
}

# Perform an M-step update for the residual covariance matrix.
update_resid_covariance_helper <- function (X, U, V, P) {
  n    <- nrow(X)
  m    <- ncol(X)
  Vnew <- matrix(0,m,m)
  for (i in 1:n) {
    out  <- compute_posterior_mvtnorm_mix(X[i,],P[i,],U,V)
    Vnew <- Vnew + out$S1 + tcrossprod(X[i,] - out$mu1)
  }
  return(Vnew/n)
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
  U <- shrink_cov(T,minval)

  # Recover the solution for the original (untransformed) data.
  return(t(R) %*% U %*% R)
}

# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance V, and z is drawn from a mixture of multivariate
# normals, each with zero mean, covariance U[,,i] and weight w[i].
# Return the posterior mean (mu1) and covariance (S1) of z. Note that
# input w1 must be the vector of *posterior* mixture weights (see
# compute_posterior_probs).
compute_posterior_mvtnorm_mix <- function (x, w1, U, V) {
  m   <- length(x)
  k   <- length(w1)
  mu1 <- rep(0,m)
  S1  <- matrix(0,m,m)
  for (i in 1:k) {
    out <- compute_posterior_mvtnorm(x,U[,,i],V)
    mu1 <- mu1 + w1[i] * out$mu1
    S1  <- S1 + w1[i] * (out$S1 + tcrossprod(out$mu1))
  }
  S1 <- S1 - tcrossprod(mu1)
  return(list(mu1 = mu1,S1 = S1))
}

# Suppose x is drawn from a multivariate normal distribution with mean
# z and covariance V, and z is drawn from a multivariate normal
# distribution with mean zero and covariance U. Return the posterior
# mean (mu) and covariance (S1) of z. These calculations will only
# work if V is positive definite (invertible).
compute_posterior_mvtnorm <- function (x, U, V) {
  m   <- length(x)
  S1  <- solve(U %*% solve(V) + diag(m)) %*% U
  mu1 <- drop(S1 %*% solve(V,x))
  return(list(mu1 = mu1,S1 = S1))
}

#' @rdname ud_fit
#'
#' @export
#' 
ud_fit_control_default <- function()
  list(weights.update = "em",   # "em" or "none"
       update.U       = "teem", # "ed", "teem" or "none"
       resid.update   = "em",   # "em" or "none"
       version        = "Rcpp", # "R" or "Rcpp"
       maxiter        = 100,
       minval         = -1e-8,
       tol            = 1e-6)
