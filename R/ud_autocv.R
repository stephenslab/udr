
# Function to perform n-fold cross-validation
#' @param X: n by R data matrix
#' @param V: R by R residual covariance matrix
#' @param nfold: An integer, the number of folds used in CV. 
#' @param n_unconstrained: An integer, the number of unconstrained matrix to fit
#' @param n_rank1: An integer, the number of rank_rank1 matrix to fit
cv_single_model = function(X, V, nfold, n_unconstrained, n_rank1, control, verbose){
  
  n = nrow(X)
  size = round(n / nfold)
  loglik.test = rep(NA, nfold)
  fit_by_fold = c()
  
    # split data
  for (i in 1:nfold){
    start = 1+(i-1)*size
    end = ifelse(i==nfold, n, i*size)
    X.test = X[start:end, ]
    X.train = X[-c(start:end), ]
    
    # fit model 
    fit1 <- ud_init(X.train, n_unconstrained = n_unconstrained, n_rank1 = n_rank1, U_scaled = NULL, V = V)
    fit2 <- ud_fit(fit1,control = list(unconstrained.update = "ted", rank_rank1.update = "ted",
                                       resid.update = 'none', maxiter = control$maxiter, tol = 1e-5),verbose = verbose)
    
    U <- lapply(fit2$U,function (e) "[["(e,"mat"))
    U <- simplify2array(U)
    loglik.test[i] <- udr:::loglik_ud(X.test, fit2$w, U, fit2$V)
    fit_by_fold[[i]] = fit2
  }
  return(list(fit_by_fold = fit_by_fold, avg.loglik = mean(loglik.test)/size))
}


# Function to select best k components by cross-validation. At least, one of k_unconstrained and k_rank1 
# should be specified by the user. 
#' @param X: n by R data matrix
#' @param V: residual covariance matrix
#' @param nfold: An integer, the number of folds used in CV. 
#' @param k_unconstrained: An integer or a vector of integers specifying the number of
#' unconstrained components to experiment with.
#' @param k_rank1: An integer or a vector of integers specifying the number of
#' rank1 components to experiment with. 

ud_fit_cv = function(X, V, nfold, k_unconstrained = 0, k_rank1= 0, control=list(), verbose){
  
  if (length(k_unconstrained) == 1 & length(k_rank1) == 1)
    if (k_unconstrained == 0 & k_rank1 == 0)
      stop("At least \"k_unconstrained\" or \"k_rank1\" has to be specified")
  
  if ((length(k_unconstrained) > 1) & (length(k_rank1) > 1))
    if (length(k_unconstrained) != length(k_rank1))
      stop("\"k_unconstrained\" and \"k_rank1\" should have the same length if they are both vectors")
  
  k = max(length(k_unconstrained), length(k_rank1))
  
  avg_logliks = c(-Inf, rep(NA, k)) # store average loglikelihood under each scenario
  fit_all = c() # store fit object
  kmat = matrix(0, nrow = 2, ncol = k) # store the values of k_unconstrained and k_rank1 under each evaluated scenario
  rownames(kmat) = c("k_unconstrained", "k_rank1")
  
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
 
  # Perform CV on different k and evaluate if average log-likelihood 
  # increases. Early stop is available if og-likelihood decreases.
  for (i in 1:k){
    n_unconstrained = ifelse(i > length(k_unconstrained), k_unconstrained, k_unconstrained[i])
    n_rank1 = ifelse(i > length(k_rank1), k_rank1, k_rank1[i])
    kmat[,i] = c(n_unconstrained, n_rank1) # store n_unconstrained and n_rank1 in curr iteration
    
    # Perform CV
    res = cv_single_model(X, V, nfold, n_unconstrained, n_rank1, control, verbose)
    # store results
    avg_logliks[i+1]= res$avg.loglik
    fit_all[[i]] = res$fit_by_fold

    diff = avg_logliks[i+1] - avg_logliks[i] # compare average loglik between curr iter and previous iter
    if (diff < 0){
      break 
    }
  }
  return(list(fit_all = fit_all, scenario = kmat[,1:i], avg_logliks= avg_logliks[2:(i+1)]))
}


#' Function to get the best_fit_cv object based on highest average test log-likelihood.
#' @param res: the result from ud_fit_cv(). 
#' @return best_fit_cv: an object storing training results, its corresponding average 
#' test loglikelihood and parameter configuration. 
get_best_fit_cv = function(res){
  
  best_fit_cv = c()
  indx = which.max(res$avg_logliks)
  best_fit_cv$res_single_model = res$fit_all[[indx]]
  best_fit_cv$test.loglik = res$avg_logliks[indx]
  best_fit_cv$n_unconstrained = unname(res$scenario[1, indx])
  best_fit_cv$n_rank1 = unname(res$scenario[2, indx])
  return(best_fit_cv)
}

#' Function to get the best fit on the whole dataset based on cv results. 
#' @param X: n by R data matrix
#' @param V: R by R residual covariance matrix
#' @param best_fit_cv: An object storing the best cv results and its parameter configuration.
#' @return best_fit: the fit object on whole dataset. 
get_best_fit = function(X, V, best_fit_cv, control = list(), verbose){
  
    control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
    # Initialize with best configuration 
    fit1 <- ud_init(X, n_unconstrained = best_fit_cv$n_unconstrained, n_rank1 = best_fit_cv$n_rank1, U_scaled = NULL, V = V)
    fit2 <- ud_fit(fit1,control = list(unconstrained.update = "ted", rank_rank1.update = "ted",
                                       resid.update = 'none', maxiter = control$maxiter, tol = 1e-5),verbose = verbose)
  return(fit2)  
}