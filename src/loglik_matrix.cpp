#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_loglik_matrix_notiid (const mat& X, const cube& U, const cube& V,
				   mat& loglik);

// FUNCTION DEFINITIONS
// --------------------
// Compute the n x k log-likelihood matrix given current estimates
// of the model parameters for the general case when the samples do
// not necessarily share the same residual covariance matrix V.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_loglik_matrix_notiid_rcpp (const arma::mat& X,
					     const arma::cube& U,
					     const arma::cube& V) {
  unsigned int n = X.n_rows;
  unsigned int k = U.n_slices;
  mat loglik(n,k);
  compute_loglik_matrix_notiid(X,U,V,loglik);
  return loglik;
}

// Compute the n x k log-likelihood matrix given current estimates of
// the model parameters for the more general case when the residual
// variance is not necessarily the same for all samples.
void compute_loglik_matrix_notiid (const mat& X, const cube& U, const cube& V, 
				   mat& loglik) {
  
  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = U.n_slices;

  // These are used to store intermediate calculations.
  mat T(m,m);
  mat L(m,m);
  vec x(m);

  // Compute the log-probabilities, stored in an n x k matrix.
  for (unsigned int j = 0; j < k; j++) {
    for (unsigned int i = 0; i < n; i++) {
      x = trans(X.row(i));
      T = V.slice(i) + U.slice(j);
      L = chol(T,"lower");
      loglik(i,j) = ldmvnorm(x,L);
    }
  }
}
