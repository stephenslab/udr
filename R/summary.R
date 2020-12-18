#' @title Summarize Ultimate Deconvolution Fit
#'
#' @description \code{summary} method for the \dQuote{ud_fit} class.
#'
#' @param object An object of class \dQuote{ud_fit}, usually the
#'   result of calling function \code{link{ud_fit}}.
#'  
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
#' 
#' @method summary ud_fit
#'
#' @return \code{summary.ud_fit} returns list of statistics
#'   summarizing the model fit.
#' 
#' @export
#'
summary.ud_fit <- function (object, ...) {
  covtypes <- sapply(object$U,function (x) attr(x,"covtype"))
  out <- list(m               = nrow(object$V),
              iid             = is.matrix(object$V),
              n_scaled        = sum(covtypes == "scaled"),
              n_rank1         = sum(covtypes == "rank1"),
              n_unconstrained = sum(covtypes == "unconstrained"),
              loglik          = object$loglik,
              numiter         = nrow(object$progress))
  class(out) <- c("summary.ud_fit","list")
  return(out)
}

#' @rdname summary.ud_fit
#'
#' @method print summary.ud_fit
#' 
#' @param x An object of class \dQuote{summary.ud_fit}, usually the
#'   result of a call to \code{summary.ud_fit}.
#'
#' @export
#' 
print.summary.ud_fit <- function (x, ...) {
  cat("Model overview:\n")
  cat(sprintf("  Dimension of data points: %d\n",x$m))
  if (x$iid)
    cat("Data points are i.i.d. (same V)\n")
  else
    cat("Data points are not i.i.d. (different Vs)\n")
  cat(sprintf("  Number of scaled prior covariances: %d\n",x$n_scaled))
  cat(sprintf("  Number of rank-1 prior covariances: %d\n",x$n_rank1))
  cat(sprintf("  Number of unconstrained prior covariances: %d\n",
              x$n_unconstrained))
  cat(sprintf("Evaluation of fit (%d updates performed):\n",x$numiter))
  cat(sprintf("  Log-likelihood: %+0.12e\n",x$loglik))
  return(invisible(x))
}
