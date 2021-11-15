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

  # Update the prior covariance matrices.
  k <- length(fit$U)
  if (is.matrix(fit$V)) {

    # Update the prior covariances in the i.i.d. case (when all the
    # residual variances V are the same).
    R <- chol(fit$V)
    Y <- fit$X %*% solve(R)
    for (i in 1:k)
      fit$U[[i]] <- do.call(covupdates[i],list(Y = Y,U = fit$U[[i]],
                                               R = R,p = fit$P[,i]))
  } else {

    # Update the prior covariances in the non-i.i.d. case (when the
    # residual variances V are not all the same).
    V <- fit$V
    if (is.list(V))
      V <- list2array(V)
    for (i in 1:k)
      fit$U[[i]] <- do.call(covupdates[i],list(X = fit$X,U = fit$U[[i]],
                                               V = V,p = fit$P[,i]))
  }
  
  # Output the updated fit.
  return(fit)
}

# Initialize the data structure for a scaled prior covariance matrix.
create_prior_covariance_struct_scaled <- function (X, U, s = 1) {
  rownames(U) <- colnames(X)
  colnames(U) <- colnames(X)
  U <- list(s = s,U0 = U,mat = U)
  attr(U,"covtype") <- "scaled"
  return(U)
}

# Initialize the data structure for a rank-1 prior covariance matrix.
create_prior_covariance_struct_rank1  <- function (X, U) {
  u <- getrank1(U)
  U <- list(vec = u,mat = tcrossprod(u))
  names(U$vec)    <- colnames(X)
  rownames(U$mat) <- colnames(X)
  colnames(U$mat) <- colnames(X)
  attr(U,"covtype") <- "rank1"
  return(U)
}

# Initialize the data structure for an unconstrained prior covariance
# matrix.
create_prior_covariance_struct_unconstrained <- function (X, U) {
  rownames(U) <- colnames(X)
  colnames(U) <- colnames(X)
  U <- list(mat = U)
  attr(U,"covtype") <- "unconstrained"
  return(U)
}

# Update the data structure for a scaled prior covariance matrix.
# Input U is the current data structure, and "s" is the new estimate
# of the scaling factor. This function is called in the
# update_prior_covariance_scaled_* functions below.
update_prior_covariance_struct_scaled <- function (U, s) {
  U$s   <- s
  U$mat <- s * U$U0
  return(U)
}

# Update the data structure for a rank-1 prior covariance matrix.
# Input U is the current data structure, and "vec" is a vector
# containing the new estimates, such that the new rank-1 matrix is
# tcrossprod(vec).
update_prior_covariance_struct_rank1 <- function (U, vec) {
  mat <- tcrossprod(vec)
  names(vec) <- names(U$vec)
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$vec <- vec
  U$mat <- mat
  return(U)
}

# Update the data structure for an unconstrained prior covariance
# matrix. Input "U" is the current data structure, and "mat" is the
# newly estimated matrix.
update_prior_covariance_struct_unconstrained <- function (U, mat) {
  rownames(mat) <- rownames(U$mat)
  colnames(mat) <- colnames(U$mat)
  U$mat <- mat
  return(U)
}

# These functions return the unconstrained prior covariance matrix
# without updating in the i.i.d. case when the residual variance V is
# the same for all data points.
update_prior_covariance_scaled_none_iid             <- function (Y,U,R,p) U
update_prior_covariance_scaled_none_iid_rcpp        <- function (Y,U,R,p) U
update_prior_covariance_rank1_none_iid              <- function (Y,U,R,p) U
update_prior_covariance_rank1_none_iid_rcpp         <- function (Y,U,R,p) U
update_prior_covariance_unconstrained_none_iid      <- function (Y,U,R,p) U
update_prior_covariance_unconstrained_none_iid_rcpp <- function (Y,U,R,p) U

# These functions return the unconstrained prior covariance matrix
# without updating it in the non-i.i.d. case when the residual
# variance V is not the same for all data points.
update_prior_covariance_scaled_none_notiid             <- function (X,U,V,p) U
update_prior_covariance_scaled_none_notiid_rcpp        <- function (X,U,V,p) U
update_prior_covariance_rank1_none_notiid              <- function (X,U,V,p) U
update_prior_covariance_rank1_none_notiid_rcpp         <- function (X,U,V,p) U
update_prior_covariance_unconstrained_none_notiid      <- function (X,U,V,p) U
update_prior_covariance_unconstrained_none_notiid_rcpp <- function (X,U,V,p) U
