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
                            ifelse(control$version == "Rcpp","_rcpp",""))
  names(covupdates) <- names(fit$U)
  return(list(control = control,covupdates = covupdates))
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

  # Get the residual covariance matrix or matrices.
  V <- fit$V
  if (is.list(V))
    V <- list2array(V)

  # Update the prior covariance matrices.
  k <- length(fit$U)
  for (i in 1:k)
    fit$U[[i]] <- do.call(covupdates[i],list(X = fit$X,U = fit$U[[i]],V = V,
                                             p = fit$P[,i],minval = minval))
  
  # Output the updated fit.
  return(fit)
}
