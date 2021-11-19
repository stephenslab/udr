#include "misc.h"
#include "posterior.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void update_resid_covariance (const mat& X, const cube& U, const mat& V,
			      const mat& P, mat& Vnew);

// FUNCTION DEFINITIONS
// --------------------
// Perform an M-step update for the residual covariance matrix.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_resid_covariance_rcpp (const arma::mat& X, 
					const arma::cube& U, 
					const arma::mat& V,
					const arma::mat& P) {
  mat Vnew = V;
  update_resid_covariance(X,U,V,P,Vnew);
  return Vnew;
}

// Perform an M-step update for the residual covariance matrix.
void update_resid_covariance (const mat& X, const cube& U, const mat& V,
			      const mat& P, mat& Vnew) {

  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;

  // These are used to store intermediate calculations.
  cube B1(m,m,k);
  mat S1(m,m);
  vec mu1(m);
  vec x(m);
  vec p(k);
  
  // Compute the posterior covariances for each mixture component. The
  // posterior covariances do not depend on x, so we compute them
  // upfront.
  for (unsigned int j = 0; j < k; j++)
    compute_posterior_covariance_mvtnorm(U.slice(j),V,B1.slice(j));

  // Compute the M-step update for the residual covariance.
  Vnew.fill(0);
  for (unsigned int i = 0; i < n; i++) {
     x = trans(X.row(i));
     p = trans(P.row(i));
     compute_posterior_mvtnorm_mix(x,p,V,B1,mu1,S1);
     mu1 -= x;
     Vnew += S1;
     Vnew += mu1 * trans(mu1);
  }
  Vnew /= n;
}
