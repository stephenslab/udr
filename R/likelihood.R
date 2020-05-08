# Compute the log-likelihood for the multivariate normal means model.
#
#' @importFrom mvtnorm dmvnorm
loglik_mvebnm <- function (X, w, U, S, version = c("Rcpp","R")) {
  version <- match.arg(version)
  if (version == "Rcpp") {
    # TO DO.
  } else if (version == "R") {
    n <- nrow(X)
    k <- length(U)
    y <- rep(0,n)
    for (i in 1:k)
      y <- y + w[i] * dmvnorm(X,sigma = S + U[[i]])
    y <- sum(log(y))
  }
  return(y)
}
