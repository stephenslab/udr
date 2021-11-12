#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void ted (const mat& X, const vec& p, mat& U, double minval, unsigned int r);

// FUNCTION DEFINITIONS
// --------------------
// Perform an M-step update for a prior covariance matrix U using the
// "eigenvalue truncation" technique described in Won et al (2013).
// Input p is a vector of weights associated with the rows of X.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ted_rcpp (const arma::mat& X, const arma::vec& p, double minval, 
		    unsigned int r) {
  unsigned int m = X.n_cols;
  mat U(m,m);
  ted(X,p,U,minval,r);
  return U;
}

// Perform an M-step update for one of the prior covariance matrices
// using the eigenvalue-truncation technique described in Won et al
// (2013).
void ted (const mat& X, const vec& p, mat& U, double minval, unsigned int r) {
  mat X1 = X;

  // Transform the data so that the residual covariance is I, then
  // compute the maximum-likelhood estimate (MLE) for T = U + I.
  vec p1 = p;
  safenormalize(p1);
  scale_rows(X1,sqrt(p1));
  mat T = crossprod(X1);

  // Find U maximizing the expected complete log-likelihood subject to
  // U being positive definite.
  shrink_cov(T,U,minval,r);
}
