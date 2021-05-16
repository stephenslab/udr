# Perform an M-step update for the residual covariance matrix (or no
# update if update = "none"). Input argument P should be the n x k
# matrix of posterior mixture assignment probabilities returned by
# compute_posterior_probs. Input argument V should be an m x m
# matrix. Input argument U may either be a list of length k in which
# U[[i]]$mat is an m x m matrix, or an m x m x k array.
#
#' @rdname ud_fit_advanced
#' 
#' @export
#' 
update_resid_covariance <- function (fit, update = c("em","none"),
                                     version = c("Rcpp","R")) {
  update  <- match.arg(update)
  version <- match.arg(version)

  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")
  if (!is.matrix(fit$V))
    stop("The residual covariance V can only be updated when it is the same ",
         "for all data points")

  # Process U.
  U <- ulist2array(fit$U)

  # Update the residual covariance matrix, V.
  if (update == "em") {
    if (version == "R")
      Vnew <- update_resid_covariance_helper(fit$X,U,fit$V,fit$P)
    else if (version == "Rcpp")
      Vnew <- update_resid_covariance_rcpp(fit$X,U,fit$V,fit$P)
  } else if (update == "none")
    Vnew <- fit$V
  else
    stop("control$resid.update == \"",update,"\" is not implemented")

  # Add back the row and column names.
  rownames(Vnew) <- rownames(fit$V)
  colnames(Vnew) <- colnames(fit$V)
  
  # Output the updated fit.
  fit$V <- Vnew
  return(fit)
}

# Implements update_resid_covariance with version = "R".
update_resid_covariance_helper <- function (X, U, V, P) {
  n <- nrow(X)
  m <- ncol(X)
  Vnew <- matrix(0,m,m)
  for (i in 1:n) {
    out <- compute_posterior_mvtnorm_mix(X[i,],P[i,],U,V)
    Vnew <- Vnew + out$S1 + tcrossprod(X[i,] - out$mu1)
  }
  return(Vnew/n)
}
