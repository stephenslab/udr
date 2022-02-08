

# Function to split data into train/test
#' @param X: n by R data matrix
data.split <- function(X){
  
  n = nrow(X)
  # split data into train/test: ratio = 4:1
  indx = sample(1:n, round(n/5), replace = FALSE)
  
  X.test = X[indx, ]
  X.train = X[-indx, ]
  return(list(X.test = X.test, X.train = X.train))
}



# Function to perform cross-validation 
#' @param X: n by R data matrix
#' @param V: residual covariance matrix
#' @param grid_k: A list of Ks to experiment with
udr_fit_cv = function(X, V, grid_k = c(1,2,3,4,6,8,10,20,30,40,50,100), control=list()){
  
  X.train = data.split(X)$X.train
  X.test = data.split(X)$X.test
  nk = length(grid_k)
  loglik.test = c(-Inf, rep(NA, n))
  
  # Train model based on different k and evaluate if 
  # test log-likelihood decreases
  for (i in 1:nk){
    
    fit1 <- ud_init(X, n_unconstrained = grid_k[i], n_rank1 = 0, U_scaled = NULL, V = V)
    fit2 <- ud_fit(fit1,control = list(unconstrained.update = "ted", rank1.update = "ted",
                                       resid.update = 'none', maxiter = 100, tol = 1e-5),verbose = FALSE)
    
    U <- lapply(fit2$U,function (e) "[["(e,"mat"))
    U <- simplify2array(U)
    loglik.test[i+1] <- udr:::loglik_ud(X.test, fit2$w, U, fit2$V)
    if (loglik.test[i+1] < loglik.test[i]){
      break
    }
  }
  return(list(k.best = grid_k[i-1], loglik.test = loglik.test))
}