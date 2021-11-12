#' @rdname ud_fit_advanced
#'
#' @title Low-Level Interface for Fitting Ultimate Deconvolution Models
#' 
#' @description These functions allow for more fine-grained fitting of
#' Ultimate Deconvolution models. Only minimal argument checking is
#' performed and this interface should only be used by experienced
#' users.
#'
#' @details \code{compute_posterior_probs} returns the matrix of
#' posterior mixture assignment probabilities ("responsibilities")
#' given current estimates of the model parameters. This implements
#' the E step in the EM algorithm.
#'
#' \code{update_mixture_weights} performs an M-step update for the
#' mixture weights in the mixture prior (or no update if \code{update
#' = "none"}).
#'
#' \code{update_resid_covariance} performs an M-step update for the
#' residual covariance matrix, V (or no update if \code{update = "none"}).
#'
#' \code{assign_prior_covariance_updates} determine the functions for
#' updating the prior covariance matrices based on the the covariance
#' matrix types and the control settings. The \code{"covupdates"}
#' return value is a character vector with one entry for each prior
#' covariance matrix.
#' 
#' \code{update_prior_covariances} performs an M-step update for all
#' the prior covariance matrices, U.
#'
#' @param fit An Ultimate Deconvolution model fit. Typically,
#'   this will be an output of \code{\link{ud_init}} or \code{ud_fit}.
#'
#' @param version When \code{version = "R"}, the R-only
#'   implementation is used; when \code{version = "Rcpp"}, the more
#'   efficient C++ implementation is used.
#' 
#' @return All functions except \code{assign_prior_covariance_updates}
#'   return an Ultimate Deconvolution model fit; see
#'   \code{\link{ud_fit}} for details.
#'
#' @seealso \code{\link{ud_init}}, \code{\link{ud_fit}}
#' 
#' @examples
#' # Simulate data from a UD model.
#' set.seed(1)
#' n <- 4000
#' V <- rbind(c(0.8,0.2),
#'            c(0.2,1.5))
#' U <- list(none   = rbind(c(0,0),
#'                          c(0,0)),
#'           shared = rbind(c(1.0,0.9),
#'                          c(0.9,1.0)),
#'           only1  = rbind(c(1,0),
#'                          c(0,0)),
#'           only2  = rbind(c(0,0),
#'                          c(0,1)))
#' w <- c(0.8,0.1,0.075,0.025)
#' rownames(V) <- c("d1","d2")
#' colnames(V) <- c("d1","d2")
#' X <- simulate_ud_data(n,w,U,V)
#' 
#' # Perform 4 EM updates.
#' fit1 <- ud_init(X)
#' fit1 <- ud_fit(fit1,control = list(maxiter = 4))
#'
#' # Update the responsibilities matrix.
#' fit2 <- compute_posterior_probs(fit1)
#'
#' # Update the mixture weights.
#' fit2 <- update_mixture_weights(fit2)
#'
#' # Update the residual covariance, V.
#' fit2 <- update_resid_covariance(fit2)
#'
#' # Update the prior covariance matrices, U, using the defaults.
#' fit2 <- update_prior_covariances(fit2)
#'
#' # Update only the unconstrained prior covariance matrices using the
#' # truncated eigenvalue decomposition ("ted") algorithm.
#' control <- list(scaled.update = "none",rank1.update = "none",
#'                 unconstrained.update = "ted")
#' updates <- assign_prior_covariance_updates(fit2,control)$covupdates
#' fit2 <- update_prior_covariances(fit2,updates)
#' 
#' # Compute the new log-likelihood and compare the old one.
#' print(logLik(fit1),digits = 8)
#' print(logLik(fit2),digits = 8)
#' 
#' @keywords internal
#' 
#' @export
#' 
compute_posterior_probs <- function (fit, version = c("Rcpp","R")) {
  version <- match.arg(version)
  
  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")
  
  # Get the prior (U) and residual (V) covariance matrices.
  U <- ulist2array(fit$U)
  V <- fit$V
  if (is.list(V))
    V <- list2array(V)

  # Compute the responsibilities matrix.
  if (is.matrix(V)) {

    # Perform the computations for the special case when the same
    # residual variance is used for all samples.
    # if (version == "R")
      fit$P <- compute_posterior_probs_iid(fit$X,fit$w,U,V)
    # else if (version == "Rcpp")
    #   fit$P <- compute_posterior_probs_iid_rcpp(fit$X,fit$w,U,V)
  } else {
      
    # Perform the computations for the more general case when the
    # residual variance is not necessarily the same for all samples.
    if (version == "R")
      fit$P <- compute_posterior_probs_general(fit$X,fit$w,U,V)
    else if (version == "Rcpp")
      fit$P <- compute_posterior_probs_general_rcpp(fit$X,fit$w,U,V)
  }

  # Add row and column names to the responsibilities matrix.
  rownames(fit$P) <- rownames(fit$X)
  colnames(fit$P) <- names(fit$U)

  # Output the updated fit.
  return(fit)
}

# This implements the calculations for compute_posterior_probs for the
# special case when the same residual covariance matrix is used for
# all samples.
#
#' @importFrom mvtnorm dmvnorm
compute_posterior_probs_iid <- function (X, w, U, V) {
      
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
compute_posterior_probs_general <- function (X, w, U, V) {
      
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

