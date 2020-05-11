#include "misc.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_mvebnm (const mat& X, const vec& w, const cube& U, 
		      const mat& S);

// FUNCTION DEFINITIONS
// --------------------
// Compute the log-likelihood for the mvebnm model. See the mvebnm R
// function for a description of the inputs.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_mvebnm_rcpp (const arma::mat& X, const arma::vec& w, 
			   const arma::cube& U, const arma::mat& S) {
  return loglik_mvebnm(X,w,U,S);
}

// Compute the log-likelihood for the mvebnm model. See the mvebnm R
// function for description of the inputs.
double loglik_mvebnm (const mat& X, const vec& w, const cube& U, 
		      const mat& S) {
  
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
    T = S + U.slice(j);
    L = chol(T,"lower");      
    for (unsigned int i = 0; i < n; i++) {
      x     = trans(X.row(i));
      y(i) += w(j) * exp(ldmvnorm(x,L));
    }
  }
  
  return(sum(log(y)));
}
