#' @rdname ud_fit_advanced
#'
#' @param control A list of parameters controlling the behaviour of
#'   the model fitting and initialization. See \code{\link{ud_fit}} for
#'   details.
#' 
#' @export
#' 
assign_prior_covariance_updates <- function (fit, control = list()) {

  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")

  # Check and process the optimization settings.
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
  if (is.na(control$scaled.update))
    control$scaled.update <- ifelse(is.matrix(fit$V),"fa","none")
  if (is.na(control$rank1.update))
    control$rank1.update  <- ifelse(is.matrix(fit$V),"ted","fa")
  if (is.na(control$unconstrained.update))
    control$unconstrained.update <- ifelse(is.matrix(fit$V),"ted","none")
    
  # Extract the "covtype" attribute from the prior covariance (U)
  # matrices.
  covtypes <- sapply(fit$U,function (x) attr(x,"covtype"))

  # Determine the names of the functions used to update the prior
  # covariance matrices.
  k <- length(covtypes)
  covupdates <- rep(as.character(NA),k)
  for (i in 1:k)
    covupdates[i] <- paste0("update_prior_covariance_",covtypes[i],"_",
                            control[[paste(covtypes[i],"update",sep = ".")]],
                            ifelse(is.matrix(fit$V),"_iid","_notiid"),
                            ifelse(control$version == "Rcpp","_rcpp",""))
  names(covupdates) <- names(fit$U)
  return(list(control = control, covupdates = covupdates))
}

#' @rdname ud_fit_advanced
#'
#' @param covupdates Functions or character strings naming the
#'   functions to be called for updating the prior covariance matrices.
#' 
#' @param minval Minimum eigenvalue allowed in the prior covariance
#'   matrices. Should be a small, positive number.
#'
#' @export
#' 
update_prior_covariances <-
  function (fit,            
            covupdates = assign_prior_covariance_updates(fit)$covupdates,
            minval = 1e-8) {

  # Check input argument "fit".
  if (!(is.list(fit) & inherits(fit,"ud_fit")))
    stop("Input argument \"fit\" should be an object of class \"ud_fit\"")

  # Update the prior covariance matrices.
  k <- length(fit$U)
  if (is.matrix(fit$V)) {
    
    # Transform the data x to x' so that x' ~ N(0,U' + I).
    fit0 <- simplify_model(fit)

    # Update the prior covariances U for the simpler model, x' ~ N(0,U' + I).
    for (i in 1:k)
      if (sum(fit0$P[,i]) == 0){
        next
      }
      fit0$U[[i]] <- do.call(covupdates[i],
                             list(X = fit0$X,U = fit0$U[[i]],
                                  p = fit0$P[,i],minval = minval))
    
    # Transform the data x' ~ N(0,U' + I) back to x ~ N(0,U + V).
    fit0 <- unsimplify_model(fit0)
    
    # Update the updated prior covariances U in original fit object
    covtypes <- sapply(fit$U,function (x) attr(x,"covtype"))
    for (i in 1:k){
      if (covtypes[i] == "scaled")
        fit$U[[i]]$s <- fit0$U[[i]]$s
      fit$U[[i]] <- fit0$U[[i]]
    }
  } else {

    # Update the prior covariances in the non-i.i.d. case (when the
    # residual variances V are not all the same). The inputs to each
    # update are: X, the data matrix; U, the current covariance matrix
    # estimate; V, the residual variances stored in an m x m x n
    # array, where n = nrow(X) and m = ncol(X); and p, a vector of
    # mixture weights associated with the rows of X. The output is the
    # new estimate of U is the model x[i,] ~ N(0,U + V[[i]]).
    V <- fit$V
    if (is.list(V))
      V <- list2array(V)
    for (i in 1:k)
      fit$U[[i]] <- do.call(covupdates[i],
                             list(X = fit$X,U = fit$U[[i]],V = V,
                                  p = fit$P[,i],minval = minval))
  }
  
  # Output the updated fit.
  return(fit)
}

# These functions return the unconstrained prior covariance matrix
# without updating in the i.i.d. case when the residual variance V is
# the same for all data points.
update_prior_covariance_scaled_none_iid             <- function (X,U,p,...) U
update_prior_covariance_scaled_none_iid_rcpp        <- function (X,U,p,...) U
update_prior_covariance_rank1_none_iid              <- function (X,U,p,...) U
update_prior_covariance_rank1_none_iid_rcpp         <- function (X,U,p,...) U
update_prior_covariance_unconstrained_none_iid      <- function (X,U,p,...) U
update_prior_covariance_unconstrained_none_iid_rcpp <- function (X,U,p,...) U

# These functions return the unconstrained prior covariance matrix
# without updating it in the non-i.i.d. case when the residual
# variances V are not the same for all data points.
update_prior_covariance_scaled_none_notiid            <-function(X,U,V,p,...) U
update_prior_covariance_scaled_none_notiid_rcpp       <-function(X,U,V,p,...) U
update_prior_covariance_rank1_none_notiid             <-function(X,U,V,p,...) U
update_prior_covariance_rank1_none_notiid_rcpp        <-function(X,U,V,p,...) U
update_prior_covariance_unconstrained_none_notiid     <-function(X,U,V,p,...) U
update_prior_covariance_unconstrained_none_notiid_rcpp<-function(X,U,V,p,...) U

# Transform data x ~ N(0,U + V) to obtain equivalent model x' ~ N(0,U' + I).
# Note that U should be a list of matrices; the transformation is
# applied to each of the matrices. The Cholesky factor R = chol(V) is
# also returned.
simplify_model <- function (fit) {
  k     <- length(fit$U)
  fit$R <- chol(fit$V)
  Rinv  <- solve(fit$R)
  fit$X <- fit$X %*% Rinv
  for (i in 1:k) {
    u <- fit$U[[i]]
    f <- paste("transform_prior_covariance_struct",attr(u,"covtype"),sep="_")
    fit$U[[i]] <- do.call(f,list(U = u,A = Rinv))
  }
  fit$V <- as.numeric(NA)
  return(fit)
}

# Undo the changes to the ud_fit object made by simplify_model(fit).
unsimplify_model <- function (fit) {
  k     <- length(fit$U)
  R     <- fit$R
  fit$X <- fit$X %*% R
  fit$V <- crossprod(fit$R)
  for (i in 1:k) {
    u <- fit$U[[i]]
    f <- paste("transform_prior_covariance_struct",attr(u,"covtype"),sep="_")
    fit$U[[i]] <- do.call(f,list(U = u,A = R))
  }
  fit$R <- as.numeric(NA)
  return(fit)
}