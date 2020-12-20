// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_posterior_probs_iid (const mat& X, const vec& w, const cube& U,
				  const mat& V, mat& P);

void compute_posterior_probs_general (const mat& X, const vec& w, 
				      const cube& U, const cube& V, 
				      mat& P);

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

void compute_posterior_mvtnorm_mix (const vec& x, const vec& w1, const mat& V,
				    const cube& S, vec& mu1, mat& S1, vec& y);

void compute_posterior_covariance_mvtnorm (const mat& U, const mat& V,
					   const mat& I, mat& S);

void compute_posterior_mean_mvtnorm (const vec& x, const mat& S, 
				     const mat& V, vec& mu1);

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
// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the special case when all the samples share the same residual
// covariance.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_posterior_probs_iid_rcpp (const arma::mat& X,
					    const arma::vec& w,
					    const arma::cube& U,
					    const arma::mat& V) {
  unsigned int n = X.n_rows;
  unsigned int k = w.n_elem;
  mat          P(n,k);
  compute_posterior_probs_iid(X,w,U,V,P);
  return P;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the more general case when the samples do not necessarily share the
// same residual covariance matrix.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_posterior_probs_general_rcpp (const arma::mat& X,
						const arma::vec& w,
						const arma::cube& U,
						const arma::cube& V) {
  unsigned int n = X.n_rows;
  unsigned int k = w.n_elem;
  mat          P(n,k);
  compute_posterior_probs_general(X,w,U,V,P);
  return P;
}

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

// Perform an M-step update for the residual covariance matrix, S.
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

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the special case when the residual variance is the same for all
// samples.
void compute_posterior_probs_iid (const mat& X, const vec& w, const cube& U,
				  const mat& V, mat& P) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;
  mat T(m,m);
  mat L(m,m);
  vec x(m);
  vec p(k);

  // Compute the log-probabilities, stored in an n x k matrix.
  for (unsigned int j = 0; j < k; j++) {
    T = V + U.slice(j);
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

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the more general case when the residual variance is not necessarily
// the same for all samples.
void compute_posterior_probs_general (const mat& X, const vec& w, 
				      const cube& U, const cube& V, 
				      mat& P) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;
  mat T(m,m);
  mat L(m,m);
  vec x(m);
  vec p(k);

  // Compute the log-probabilities, stored in an n x k matrix.
  for (unsigned int j = 0; j < k; j++) {
    for (unsigned int i = 0; i < n; i++) {
      x      = trans(X.row(i));
      T      = V.slice(i) + U.slice(j);
      L      = chol(T,"lower");
      P(i,j) = log(w(j)) + ldmvnorm(x,L);
    }
  }

  // Normalize the probabilities so that each row of P sums to 1.
  for (unsigned int i = 0; i < n; i++)
    P.row(i) = softmax(P.row(i));
}

// Perform an M-step update for the prior covariance matrices using the
// update formla derived in Bovy et al (2011).
void update_prior_covariances_ed (const mat& X, cube& U, const mat& V, 
				  const mat& P) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  mat X1(n,m);
  mat T(m,m);
  mat B(m,m);
  for (unsigned int i = 0; i < k; i++) {
    X1 = X;
    update_prior_covariance_ed(X1,U.slice(i),V,P.col(i),T,B);
  }
}

// Perform an M-step update for the prior covariance matrices, using
// the eigenvalue-truncation technique described in Won et al (2013).
void update_prior_covariances_teem (const mat& X, const mat& V, const mat& P, 
				    cube& U, double minval) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;
  mat R = chol(V,"upper");
  mat X1(n,m);
  mat T(m,m);
  mat Y(m,m);
  vec d(m);
  for (unsigned int i = 0; i < k; i++) {
    X1 = X;
    update_prior_covariance_teem(X1,R,P.col(i),U.slice(i),T,Y,d,minval);
  }
}

// Perform an M-step update for one of the prior covariance matrices
// using the update formula derived in Bovy et al (2011). 
// 
// Note that data matrix X is modified to perform the update, so
// should not be reused.
void update_prior_covariance_ed (mat& X, mat& U, const mat& V, 
				 const vec& p, mat& T, mat& B) {
  scale_rows(X,sqrt(p/sum(p)));
  T  = V + U;
  B  = solve(T,U);
  X *= B;
  U += crossprod(X) - U*B;
}

// Perform an M-step update for one of the prior covariance matrices
// using the eigenvalue-truncation technique described in Won et al
// (2013). 
//
// Input R should be R = chol(V,"upper"). Inputs T and Y are matrices
// of the same dimension as U and V storing intermediate calculations,
// and d is a vector of length m also storing an intermediate result.
//
// Note that data matrix X is modified to perform the update, so
// should not be reused.
void update_prior_covariance_teem (mat& X, const mat& R, const vec& p, 
				   mat& U, mat& T, mat& Y, vec& d, 
				   double minval) {

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
  shrink_cov(T,U,Y,d,minval);

  // Recover the solution for the original (untransformed) data.
  U = trans(R) * U * R;
}

// Perform an M-step update for the residual covariance matrix, V.
void update_resid_covariance (const mat& X, const cube& U, const mat& V,
			      const mat& P, mat& Vnew) {
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = P.n_cols;

  // Compute the posterior covariances for each mixture component. The
  // posterior covariances do not depend on x, so we compute them
  // upfront.
  cube S1(m,m,k);
  mat  I(m,m,fill::eye);
  for (unsigned int j = 0; j < k; j++)
    compute_posterior_covariance_mvtnorm(U.slice(j),V,I,S1.slice(j));

  // Compute the M-step update for the residual covariance.
  vec x(m);
  vec mu1(m);
  vec y(m);
  vec p(k);
  mat S(m,m);
  Vnew.fill(0);
  for (unsigned int i = 0; i < n; i++) {
     x = trans(X.row(i));
     p = trans(P.row(i));
     compute_posterior_mvtnorm_mix(x,p,V,S1,mu1,S,y);
     mu1  -= x;
     Vnew += S1;
     Vnew += mu1 * trans(mu1);
  }
  Vnew /= n;
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a mixture of
// multivariate normals, each with zero mean, covariance U[,,i] and
// weight w[i]. Return the posterior mean (mu1) and covariance (S1) of
// z. Note that input w1 must be the vector of *posterior* mixture
// weights (see compute_posterior_probs), and S[,,i] should be the
// posterior covariance matrix for mixture component i. Input y is
// used to store intermediate calculations; it is a vector of the same
// size as x.
void compute_posterior_mvtnorm_mix (const vec& x, const vec& w1, const mat& V,
				    const cube& S, vec& mu1, mat& S1, vec& y) {
  unsigned int k = w1.n_elem;
  mu1.fill(0);
  S1.fill(0);
  for (unsigned int i = 0; i < k; i++) {
    compute_posterior_mean_mvtnorm(x,S.slice(i),V,y);
    mu1 += w1(i) * y;
    S1  += w1(i) * (S.slice(i) + y * trans(y));
  }
  S1 -= mu1 * trans(mu1);
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a multivariate normal
// distribution with mean zero and covariance U. Return in S the
// posterior covariance of z. (Note that the posterior covariance does
// not depend on x, so it is not one of the inputs.) These calculations
// will only work if V is positive definite (invertible).
//
// Input I should be the identity matrix of the same dimension as U
// and V.
//
void compute_posterior_covariance_mvtnorm (const mat& U, const mat& V,
					   const mat& I, mat& S) {
  S = inv(U*inv(V) + I)*U;
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a multivariate normal
// distribution with mean zero and covariance U. Given S, the previously
// calculated posterior covariance of z, return the posterior mean
// (mu1) of z.
//
// Input I should be the identity matrix of the same dimension as U
// and V.
//
void compute_posterior_mean_mvtnorm (const vec& x, const mat& S, 
				     const mat& V, vec& mu1) {
  mu1 = S * solve(V,x);
}

// "Shrink" matrix T = U + I; that is, find the "best" matrix T
// satisfying the constraint that U is positive definite. This is
// achieved by setting any eigenvalues of T less than 1 to 1 + minval,
// or, equivalently, setting any eigenvalues of U less than 0 to be
// minval.
//
// Inputs d and Y are used to store the eigenvalue decomposition of T;
// these will store the eigenvalues and eigenvectors, respectively.
void shrink_cov (const mat& T, mat& U, mat& Y, vec& d, double minval) {
  unsigned int m = T.n_rows;
  if (m == 1)
    
    // Handle univariate case (m = 1).
    U(0) = max(T(0) - 1,minval);
  else {
    eig_sym(d,Y,T);
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
}
