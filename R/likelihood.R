#' @title Ultimate Deconvolution Model Likelihoods
#' 
#' @description Compute the log-likelihood for the Ultimate
#' Deconvolution model.
#'
#' @param object An Ultimate Deconvolution model fit. Typically,
#'   this will be an output of \code{\link{ud_init}} or \code{ud_fit}.
#'
#' @param version When \code{version == "R"}, the computations are
#'   performed entirely in R; when \code{version == "Rcpp"}, an Rcpp
#'   implementation is used.
#'
#' @param \dots Additional arguments (unused).
#' 
#' @return A number giving the log-likelihood for the model.
#'
#' @seealso \code{\link{ud_init}}, \code{\link{ud_fit}}
#' 
#' @method logLik ud_fit
#' 
#' @importFrom stats logLik
#'
#' @export
#' 
logLik.ud_fit <- function (object, ...) {
  if (!(is.list(object) & inherits(object,"ud_fit")))
    stop("Input argument \"object\" should be an object of class \"ud_fit\"")
  out <- sum(apply(t(object$loglik_matrix)+log(object$w),2,log_sum_exp))
  class(out) <- "logLik"
  attr(out,"df") <- as.numeric(NA)
  return(out)
}

compute_loglik <- function(fit){
  fit$loglik = as.numeric(logLik(fit))
  return(fit)
}