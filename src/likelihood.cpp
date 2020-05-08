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

  // int n = X_mat.n_rows;
  // int k = w_vec.size();

  // arma::vec y = arma::zeros<arma::vec>(n);
  // for (unsigned j = 0; j < k; ++j) {
  //   y = y + w_vec(j) * dmvnorm_mat(trans(X_mat), 
  // 				   arma::zeros<arma::vec>(X_mat.n_cols),
  // 				   T_cube.slice(j));
  // }

  // return(sum(log(y)));
  return(0);
}
