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

// TO DO: Explain here what this function does, and how to use it.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_mvebnm_rcpp (const arma::mat& X, const arma::vec& w, 
			   const arma::cube& U, const arma::mat& S) {
  return loglik_mvebnm(X,w,U,S);
}

// Compute the log-probability density of the multivariate normal
// distribution with zero mean and covariance matrix S.
double ldmvnorm (const vec& x, const mat& S) {
  mat    L = chol(S,"lower");
  double d = norm(solve(L,x),2);
  return -d*d/2 - sum(log(sqrt(2*M_PI)*L.diag()));
}

// TO DO: Explain here what this function does, and how to use it.
double loglik_mvebnm (const mat& X, const vec& w, const cube& U, 
		      const mat& S) {
  
  // Get the number of rows (n) and columns (m) of X, and the number
  // of mixture components (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;

  double y = 0;
  double t;
  mat    T(m,m);
  vec    x(m);

  // Compute the log-likelihood for each sample (row of X).
  for (unsigned int i = 0; i < n; i++) {
    x = trans(X.row(i));
    t = 0;
    for (unsigned int j = 0; j < k; j++) {
      T  = S + U.slice(j);
      t += w(j) * exp(ldmvnorm(x,T));
    }
    y += log(t);
  }
  return(y);
}
