#' @title Summarize Ultimate Deconvolution Fit
#'
#' @description Description goes here.
#'
#' @param object Describe input argument "object" here.
#'  
#' @param \dots Additional arguments passed to the generic \code{summary}
#'   or \code{print.summary} method.
#' 
#' @method summary ud_fit
#' 
#' @export
#'
summary.ud_fit <- function (object, ...) {
  out <- list()
  # TO DO.
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
  cat("A ud_fit object.\n")
  return(invisible(x))
}
