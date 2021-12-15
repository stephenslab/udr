#include "misc.h"

using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
double loglik_ud_iid (const mat& X, const vec& w, const cube& U, const mat& V);

double loglik_ud_notiid (const mat& X, const vec& w, const cube& U, 
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
// when the residual covariance V is not necessarily the same for all
// samples.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double loglik_ud_notiid_rcpp (const arma::mat& X, const arma::vec& w, 
			       const arma::cube& U, const arma::cube& V) {
  return loglik_ud_notiid(X,w,U,V);
}

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance is the same for all samples.
double loglik_ud_iid (const mat& X, const vec& w, const cube& U, 
		      const mat& V) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  vec y(n);
  vec x(m);
  for (unsigned int i = 0; i < n; i++) {
    x = trans(X.row(i));
    y(i) = ldmvnormmix(x,w,U,V);
  }
  
  return(sum(y));
}

// Compute the log-likelihood for the Ultimate Deconvolution model
// when the residual covariance is not necessarily the same for all
// samples.
double loglik_ud_notiid (const mat& X, const vec& w, const cube& U, 
			 const cube& V) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  vec y(n);
  vec x(m);
  for (unsigned int i = 0; i < n; i++) {
    x = trans(X.row(i));
    y(i) = ldmvnormmix(x,w,U,V.slice(i));
  }
  return(sum(y));
}
