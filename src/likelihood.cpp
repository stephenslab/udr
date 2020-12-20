#include "misc.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_ud_iid (const mat& X, const vec& w, const cube& U, 
		      const mat& V);

double loglik_ud_general (const mat& X, const vec& w, const cube& U, 
                          const cube& V);

// FUNCTION DEFINITIONS
// --------------------

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance V is the same for all samples.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_ud_iid_rcpp (const arma::mat& X, const arma::vec& w, 
			   const arma::cube& U, const arma::mat& V) {
  return loglik_ud_iid(X,w,U,V);
}

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance V is *not* the same for all samples.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_ud_general_rcpp (const arma::mat& X, const arma::vec& w, 
			       const arma::cube& U, const arma::cube& V) {
  return loglik_ud_general(X,w,U,V);
}

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance is the same for all samples.
double loglik_ud_iid (const mat& X, const vec& w, const cube& U, 
		      const mat& V) {
  
  // Get the number of rows (n) and columns (m) of X, and the number
  // of mixture components (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;

  mat T(m,m);
  mat L(m,m);
  vec x(m);
  vec y(n,fill::zeros);
  
  // Compute the log-likelihood for each sample (row of X).
  for (unsigned int j = 0; j < k; j++) {
    T = V + U.slice(j);
    L = chol(T,"lower");      
    for (unsigned int i = 0; i < n; i++) {
      x     = trans(X.row(i));
      y(i) += w(j) * exp(ldmvnorm(x,L));
    }
  }
  
  return(sum(log(y)));
}

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance is *not* the same for all samples.
double loglik_ud_general (const mat& X, const vec& w, const cube& U, 
			  const cube& V) {
  
  // Get the number of rows (n) and columns (m) of X, and the number
  // of mixture components (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;

  mat T(m,m);
  mat L(m,m);
  vec x(m);
  vec y(n,fill::zeros);
  
  // Compute the log-likelihood for each sample (row of X).
  for (unsigned int j = 0; j < k; j++)
    for (unsigned int i = 0; i < n; i++) {
      x     = trans(X.row(i));
      T     = V.slice(i) + U.slice(j);
      L     = chol(T,"lower");      
      y(i) += w(j) * exp(ldmvnorm(x,L));
    }
  
  return(sum(log(y)));
}
