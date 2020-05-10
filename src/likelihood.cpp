#include "likelihood.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// This is mainly used for testing the ldmvnorm C++ function.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ldmvnorm_rcpp (const arma::vec& x, const arma::mat& S) {
  return ldmvnorm(x,S);
}

// Compute the log-likelihood for the mvebnm model. See the mvebnm R
// function for a description of the inputs.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_mvebnm_rcpp (const arma::mat& X, const arma::vec& w, 
			   const arma::cube& U, const arma::mat& S) {
  return loglik_mvebnm(X,w,U,S);
}

// Compute the log-probability of x, where x is multivariate normal
// with mean zero and covariance matrix S. Input argument should be L
// be the Cholesky factor of S; L = chol(S,"lower").
double ldmvnorm (const vec& x, const mat& L) {
  double d = norm(solve(L,x),2);
  return -d*d/2 - sum(log(sqrt(2*M_PI)*L.diag()));
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
