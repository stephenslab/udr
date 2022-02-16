
# Function to perform n-fold cross-validation
#' @param X: n by R data matrix
#' @param V: R by R residual covariance matrix
#' @param nfold: An integer, the number of folds used in CV. 
#' @param n_unconstrained: An integer, the number of unconstrained matrix to fit
#' @param n_rank1: An integer, the number of rank1 matrix to fit
ud_cv = function(X, V, nfold, n_unconstrained, n_rank1, control, verbose){
  
  n = nrow(X)
  size = round(n / nfold)
  loglik.test = rep(NA, nfold)
  
  for (i in 1:nfold){
    # split data
    start = 1+(i-1)*size
    end = ifelse(i==nfold, n, i*size)
    X.test = X[start:end, ]
    X.train = X[-c(start:end), ]
    
    # fit model 
    fit1 <- ud_init(X.train, n_unconstrained = n_unconstrained, n_rank1 = n_rank1, U_scaled = NULL, V = V)
    fit2 <- ud_fit(fit1,control = list(unconstrained.update = "ted", rank1.update = "ted",
                                       resid.update = 'none', maxiter = control$maxiter, tol = 1e-5),verbose = verbose)
    
    U <- lapply(fit2$U,function (e) "[["(e,"mat"))
    U <- simplify2array(U)
    loglik.test[i] <- udr:::loglik_ud(X.test, fit2$w, U, fit2$V)
  }
  return(avg.loglik = mean(loglik.test)/size)
}


# Function to select best k components by cross-validation. At least, one of ku and k1 
# should be specified by the user. 
#' @param X: n by R data matrix
#' @param V: residual covariance matrix
#' @param nfold: An integer, the number of folds used in CV. 
#' @param ku: An integer or a list of integers specifying the number of
#' unconstrained components to experiment with.
#' @param k1: An integer or a list of integers specifying the number of
#' rank1 components to experiment with. 

ud_fit_cv = function(X, V, nfold, ku = 0, k1= 0, control=list(), verbose){
  
  k = max(length(ku), length(k1))
  avg_logliks = c(-Inf, rep(NA, k)) # store average loglikelihood under each scenario
  kmat = matrix(0, nrow = 2, ncol = k) # store the values of ku and k1 under each evaluated scenario
  rownames(kmat) = c("ku", "k1")
  
  control <- modifyList(ud_fit_control_default(),control,keep.null = TRUE)
 
  # Perform CV on different k and evaluate if average log-likelihood 
  # increases. Early stop is available if og-likelihood decreases. 
  for (i in 1:k){
    n_unconstrained = ifelse(i > length(ku), rev(ku)[1], ku[i])
    n_rank1 = ifelse(i > length(k1), rev(k1)[1], k1[i])
    kmat[,i] = c(n_unconstrained, n_rank1) # store n_unconstrained and n_rank1 in curr iteration
    
    # Perform CV
    avg_logliks[i+1]= ud_cv(X, V, nfold, n_unconstrained, n_rank1, control, verbose)
    diff = avg_logliks[i+1] - avg_logliks[i] # compare average loglik between curr iter and previous iter
    
    if (diff < 0){
      n_unconstrained = kmat[1, i-1]  # obtain n_unconstrained/n_rank1 from prev iter
      n_rank1 = kmat[2, i-1]
      break 
    }
  }
  return(list(ku = n_unconstrained, k1 = n_rank1, scenario = kmat[,1:i], avg_logliks= avg_logliks[2:i]))
}