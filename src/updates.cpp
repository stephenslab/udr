#include "misc.h"
#include "posterior.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void update_prior_covariances_ed (const mat& X, cube& U, const mat& V, 
				  const mat& P);

void update_prior_covariances_teem (const mat& X, const mat& V, const mat& P, 
				    cube& U, double minval);

void update_prior_covariance_ed (mat& X, mat& U, const mat& V, 
				 const vec& p, mat& T, mat& B);

void update_prior_covariance_teem (mat& X, const mat& R, const vec& p, mat& U, 
				   mat& T, mat& Y, vec& d, double minval);

void update_resid_covariance (const mat& X, const cube& U, const mat& V,
			      const mat& P, mat& Vnew);

void shrink_cov (const mat& T, mat& U, mat& Y, vec& d, double minval);

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Return a or b, which ever is larger.
inline double max (double a, double b) {
  double y;
  if (a > b)
    y = a;
  else
    y = b;
  return y;
}

// FUNCTION DEFINITIONS
// --------------------
// Perform an M-step update for the prior covariance matrices using the
// update formula derived in Bovy et al (2011).
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_prior_covariances_ed_rcpp (const arma::mat& X, 
					     const arma::cube& U, 
					     const arma::mat& V, 
					     const arma::mat& P) {
  cube Unew = U;
  update_prior_covariances_ed(X,Unew,V,P);
  return Unew;
}

// Perform an M-step update for the prior covariance matrices using
// the eigenvalue-truncation technique described in Won et al (2013).
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_prior_covariances_teem_rcpp (const arma::mat& X, 
					       const arma::mat& V,
					       const arma::mat& P,
					       double minval) {
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  cube U(m,m,k);
  update_prior_covariances_teem(X,V,P,U,minval);
  return U;
}

// Perform an M-step update for the residual covariance matrix.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_resid_covariance_rcpp (const arma::mat& X, 
					const arma::cube& U, 
					const arma::mat& V,
					const arma::mat& P) {
  unsigned int m = X.n_cols;
  mat Vnew(m,m);
  update_resid_covariance(X,U,V,P,Vnew);
  return Vnew;
}

// Perform an M-step update for the (unconstrained) prior covariance
// matrices using the update formula derived in Bovy et al (2011).
void update_prior_covariances_ed (const mat& X, cube& U, const mat& V, 
				  const mat& P) {

  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;

  // These matricecs are used to store intermediate calculations.
  mat X1(n,m);
  mat T(m,m);
  mat B(m,m);

  // Repeat for each prior covariance matrix to update.
  for (unsigned int i = 0; i < k; i++) {
    X1 = X;
    update_prior_covariance_ed(X1,U.slice(i),V,P.col(i),T,B);
  }
}

// Perform an M-step update for the (unconstrained) prior covariance
// matrices, using the eigenvalue-truncation technique described in
// Won et al (2013).
void update_prior_covariances_teem (const mat& X, const mat& V, const mat& P, 
				    cube& U, double minval) {

  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;

  // These are used to store intermediate calculations.
  mat R = chol(V,"upper");
  mat X1(n,m);
  mat T(m,m);
  mat Y(m,m);
  vec d(m);

  // Repeat for each prior covariance matrix to update.
  for (unsigned int i = 0; i < k; i++) {
    X1 = X;
    update_prior_covariance_teem(X1,R,P.col(i),U.slice(i),T,Y,d,minval);
  }
}

// Perform an M-step update for one of the prior covariance matrices
// using the update formula derived in Bovy et al (2011). Note that
// data matrix X is modified to perform the update, so should not be
// reused.
void update_prior_covariance_ed (mat& X, mat& U, const mat& V, const vec& p,
				 mat& T, mat& B) {
  vec p1 = p;
  safenormalize(p1);
  scale_rows(X,sqrt(p1));
  T = V + U;
  B = solve(T,U);
  X *= B;
  U += crossprod(X) - U*B;
}

// Perform an M-step update for one of the prior covariance matrices
// using the eigenvalue-truncation technique described in Won et al
// (2013). Input R should be R = chol(V,"upper"). Inputs T and Y are
// matrices of the same dimension as U and V storing intermediate
// calculations, and d is a vector of length m also storing an
// intermediate result. Also note that data matrix X is modified to
// perform the update, so should not be reused.
void update_prior_covariance_teem (mat& X, const mat& R, const vec& p, 
				   mat& U, mat& T, mat& Y, vec& d, 
				   double minval) {

  // Transform the data so that the residual covariance is I, then
  // compute the maximum-likelhood estimate (MLE) for T = U + I.
  vec p1 = p;
  safenormalize(p1);
  scale_rows(X,sqrt(p1));
  X *= inv(R);
  T  = crossprod(X);

  // Find U maximizing the expected complete log-likelihood subject to
  // U being positive definite. This update for U is based on the fact
  // that the covariance matrix that minimizes the likelihood subject
  // to the constraint that U is positive definite is obtained by
  // truncating the eigenvalues of T = U + I less than 1 to be 1; see
  // Won et al (2013), p. 434, the sentence just after equation (16).
  shrink_cov(T,U,Y,d,minval);

  // Recover the solution for the original (untransformed) data.
  U = trans(R) * U * R;
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
  mat I(m,m,fill::eye);
  mat S1(m,m);
  vec mu1(m);
  vec x(m);
  vec y(m);
  vec p(k);
  
  // Compute the posterior covariances for each mixture component. The
  // posterior covariances do not depend on x, so we compute them
  // upfront.
  for (unsigned int j = 0; j < k; j++)
    compute_posterior_covariance_mvtnorm(U.slice(j),V,I,B1.slice(j));

  // Compute the M-step update for the residual covariance.
  Vnew.fill(0);
  for (unsigned int i = 0; i < n; i++) {
     x = trans(X.row(i));
     p = trans(P.row(i));
     compute_posterior_mvtnorm_mix(x,p,V,B1,mu1,S1,y);
     mu1 -= x;
     Vnew += S1;
     Vnew += mu1 * trans(mu1);
  }
  Vnew /= n;
}

// "Shrink" matrix T = U + I; that is, find the "best" matrix T
// satisfying the constraint that U is positive definite. This is
// achieved by setting any eigenvalues of T less than 1 to 1 + minval,
// or, equivalently, setting any eigenvalues of U less than 0 to be
// minval. Inputs d and Y are used to store the eigenvalue
// decomposition of T; these will store the eigenvalues and
// eigenvectors, respectively.
void shrink_cov (const mat& T, mat& U, mat& Y, vec& d, double minval) {
  unsigned int m = T.n_rows;
  eig_sym(d,Y,T);
  minval = max(0,minval);
  for (unsigned int i = 0; i < m; i++)
    d(i) = max(d(i) - 1,minval);

  // These next few lines are equivalent to
  //
  //   U = Y * diagmat(d) * trans(Y)
  //
  // but implemented in a slightly more efficient way because they
  // avoid an extra matrix-matrix multiplication.
  inplace_trans(Y);
  scale_rows(Y,sqrt(d));
  U = crossprod(Y);
}
