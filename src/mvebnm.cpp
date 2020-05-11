// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_posterior_probs (const mat& X, const vec& w, const cube& U,
			      const mat& S, mat& P);

void update_prior_covariances_ed (const mat& X, cube& U, const mat& S, 
				  const mat& P);

void update_prior_covariances_teem (const mat& X, const mat& S, const mat& P, 
				    cube& U, double e);

void update_prior_covariance_ed (mat& X, mat& U, const mat& S, 
				 const vec& p, mat& T, mat& B);

void update_prior_covariance_teem (mat& X, const mat& R, const vec& p, mat& U, 
				   mat& T, mat& V, vec& d, double e);

void shrink_cov (const mat& T, mat& U, mat& V, vec& d, double e);

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
// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_posterior_probs_rcpp (const arma::mat& X,
					const arma::vec& w,
					const arma::cube& U,
					const arma::mat& S) {
  unsigned int n = X.n_rows;
  unsigned int k = w.n_elem;
  mat          P(n,k);
  compute_posterior_probs(X,w,U,S,P);
  return P;
}

// Perform an M-step update for the prior covariance matrices using the
// update formla derived in Bovy et al (2011).
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_prior_covariances_ed_rcpp (const arma::mat& X, 
					     const arma::cube& U, 
					     const arma::mat& S, 
					     const arma::mat& P) {
  cube Unew = U;
  update_prior_covariances_ed(X,Unew,S,P);
  return Unew;
}

// Perform an M-step update for the covariance matrices in the
// mixture-of-multivariate-normals prior.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube update_prior_covariances_teem_rcpp (const arma::mat& X, 
					       const arma::mat& S,
					       const arma::mat& P,
					       double e) {
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  cube U(m,m,k);
  update_prior_covariances_teem(X,S,P,U,e);
  return U;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters.
void compute_posterior_probs (const mat& X, const vec& w, const cube& U,
			      const mat& S, mat& P) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;
  mat T(m,m);
  mat L(m,m);
  vec x(m);
  vec p(k);

  // Compute the log-probabilities, stored in an n x k matrix.
  for (unsigned int j = 0; j < k; j++) {
    T = S + U.slice(j);
    L = chol(T,"lower");
    for (unsigned int i = 0; i < n; i++) {
      x      = trans(X.row(i));
      P(i,j) = log(w(j)) + ldmvnorm(x,L);
    }
  }

  // Normalize the probabilities so that each row of P sums to 1.
  for (unsigned int i = 0; i < n; i++)
    P.row(i) = softmax(P.row(i));
}

// Perform an M-step update for the prior covariance matrices using the
// update formla derived in Bovy et al (2011).
void update_prior_covariances_ed (const mat& X, cube& U, const mat& S, 
				  const mat& P) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  mat Y(n,m);
  mat T(m,m);
  mat B(m,m);
  for (unsigned int i = 0; i < k; i++) {
    Y = X;
    update_prior_covariance_ed(Y,U.slice(i),S,P.col(i),T,B);
  }
}

// Perform an M-step update for the prior covariance matrices, using
// the eigenvalue-truncation technique described in Won et al (2013).
void update_prior_covariances_teem (const mat& X, const mat& S, const mat& P, 
				    cube& U, double e) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  mat R = chol(S,"upper");
  mat Y(n,m);
  mat T(m,m);
  mat V(m,m);
  vec d(m);
  for (unsigned int i = 0; i < k; i++) {
    Y = X;
    update_prior_covariance_teem(Y,R,P.col(i),U.slice(i),T,V,d,e);
  }
}

// Perform an M-step update for one of the prior covariance matrices
// using the update formula derived in Bovy et al (2011). 
void update_prior_covariance_ed (mat& X, mat& U, const mat& S, 
				 const vec& p, mat& T, mat& B) {
  scale_rows(X,sqrt(p/sum(p)));
  T  = S + U;
  B  = solve(T,U);
  X *= B;
  U += crossprod(X) - U*B;
}

// For testing only.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_prior_covariance_teem_rcpp (const arma::mat& X, 
					     const arma::mat& S,
					     const arma::vec& p,
					     double e) {
  unsigned int m = X.n_cols;
  mat U(m,m);
  mat T(m,m);
  mat V(m,m);
  vec d(m);
  mat Y = X;
  mat R = chol(S,"upper");
  update_prior_covariance_teem(Y,R,p,U,T,V,d,e);
  return U;
}

// Perform an M-step update for one of the prior covariance matrices
// using the eigenvalue-truncation technique described in Won et al
// (2013). Input R should be R = chol(S,"upper"). 
void update_prior_covariance_teem (mat& X, const mat& R, const vec& p, 
  mat& U, mat& T, mat& V, vec& d, double e) {

  // Transform the data so that the residual covariance is I, then
  // compute the maximum-likelhood estimate (MLE) for T = U + I.
  scale_rows(X,sqrt(p/sum(p)));
  X *= inv(R);
  T  = crossprod(X);

  // Find U maximizing the expected complete log-likelihood subject to
  // U being positive definite. This update for U is based on the fact
  // that the covariance matrix that minimizes the likelihood subject
  // to the constraint that U is positive definite is obtained by
  // truncating the eigenvalues of T = U + I less than 1 to be 1; see
  // Won et al (2013), p. 434, the sentence just after equation (16).
  shrink_cov(T,U,V,d,e);

  // Recover the solution for the original (untransformed) data.
  U = trans(R) * U * R;
}

// "Shrink" matrix T = U + I; that is, find the "best" matrix T
// satisfying the constraint that U is positive definite. This is
// achieved by setting any eigenvalues of T less than 1 to 1 + e,
// or, equivalently, setting any eigenvalues of U less than 0 to be e.
void shrink_cov (const mat& T, mat& U, mat& V, vec& d, double e) {
  unsigned int m = T.n_rows;
  eig_sym(d,V,T);
  for (unsigned int i = 0; i < m; i++)
    d(i) = max(d(i) - 1,e);
  U = V * diagmat(d) * trans(V);
}
