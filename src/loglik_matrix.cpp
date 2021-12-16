#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_loglik_matrix_iid (const mat& X, const cube& U,
				  const mat& V, mat& P);

void compute_loglik_matrix_notiid (const mat& X,
				     const cube& U, const cube& V, 
				     mat& P);

// FUNCTION DEFINITIONS
// --------------------
// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the special case when all the samples share the same residual
// covariance V.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_loglik_matrix_iid_rcpp (const arma::mat& X,
					    const arma::cube& U,
					    const arma::mat& V) {
  unsigned int n = X.n_rows;
  unsigned int k = U.n_slices;
  mat P(n,k);
  compute_loglik_matrix_iid(X,U,V,P);
  return P;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the more general case when the samples do not necessarily share the
// same residual covariance matrix V.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_loglik_matrix_notiid_rcpp (const arma::mat& X,
					       const arma::cube& U,
					       const arma::cube& V) {
  unsigned int n = X.n_rows;
  unsigned int k = U.n_slices;
  mat P(n,k);
  compute_loglik_matrix_notiid(X,U,V,P);
  return P;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the special case when the residual variance is the same for all
// samples.
void compute_loglik_matrix_iid (const mat& X, const cube& U,
				  const mat& V, mat& P) {
  
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
    T = V + U.slice(j);
    L = chol(T,"lower");
    for (unsigned int i = 0; i < n; i++) {
      x = trans(X.row(i));
      P(i,j) = ldmvnorm(x,L);
    }
  }

}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the more general case when the residual variance is not necessarily
// the same for all samples.
void compute_loglik_matrix_notiid (const mat& X, const cube& U, const cube& V, 
				     mat& P) {
  
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
      P(i,j) = ldmvnorm(x,L);
    }
  }

}

