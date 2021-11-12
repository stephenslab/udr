#include "misc.h"
#include "posterior.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void update_prior_covariance_ed_iid (const mat& X, mat& U, const mat& V, 
				     const vec& p);

void update_resid_covariance (const mat& X, const cube& U, const mat& V,
			      const mat& P, mat& Vnew);

// FUNCTION DEFINITIONS
// --------------------
// Perform an M-step update for a prior covariance matrix U using
// the update formula derived in Bovy et al (2011), for the special
// case when the residual covariances V are the same for all data
// points.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_prior_covariance_ed_iid_rcpp (const arma::mat& X, 
					       const arma::mat& U, 
					       const arma::mat& V, 
					       const arma::vec& p) {
  mat Unew = U;
  update_prior_covariance_ed_iid(X,Unew,V,p);
  return Unew;
}

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

// Perform an M-step update for one of the prior covariance matrices
// using the update formula derived in Bovy et al (2011), for the
// special case when all the residual covariance matrices V are all
// the same.
void update_prior_covariance_ed_iid (const mat& X, mat& U, const mat& V, 
				     const vec& p) {
  vec p1 = p;
  mat X1 = X;
  safenormalize(p1);
  scale_rows(X1,sqrt(p1));
  mat T = V + U;
  mat B = solve(T,U);
  X1 *= B;
  U += crossprod(X1) - U*B;
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

