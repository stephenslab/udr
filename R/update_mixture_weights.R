#' @rdname ud_fit_advanced
#' 
#' @param update The method to use for updating the prior mixture
#'  weights or the residual covariance matrix; see \code{\link{ud_fit}}
#'  for details.
#' 
#' @export
#' 
update_mixture_weights <- function (fit, update = c("em","none")) {
  update <- match.arg(update)
    
  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")
    
  # Update the mixture weights.
  if (update == "em")
    wnew <- colSums(fit$P)/nrow(fit$P)
  else if (update == "none")
    wnew <- fit$w
  else 
    stop("control$weights.update == \"",update,"\" is not implemented")

  # Add the names back.
  names(wnew) <- names(fit$w)
  
  # Output the updated fit.
  fit$w <- wnew
  return(fit)
}

