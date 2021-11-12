#include "misc.h"
#include "posterior.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void compute_posterior_probs_iid (const mat& X, const vec& w, const cube& U,
				  const mat& V, mat& P);

void compute_posterior_probs_notiid (const mat& X, const vec& w, 
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
arma::mat compute_posterior_probs_iid_rcpp (const arma::mat& X,
					    const arma::vec& w,
					    const arma::cube& U,
					    const arma::mat& V) {
  unsigned int n = X.n_rows;
  unsigned int k = w.n_elem;
  mat P(n,k);
  compute_posterior_probs_iid(X,w,U,V,P);
  return P;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the more general case when the samples do not necessarily share the
// same residual covariance matrix V.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat compute_posterior_probs_notiid_rcpp (const arma::mat& X,
					       const arma::vec& w,
					       const arma::cube& U,
					       const arma::cube& V) {
  unsigned int n = X.n_rows;
  unsigned int k = w.n_elem;
  mat P(n,k);
  compute_posterior_probs_notiid(X,w,U,V,P);
  return P;
}

// Compute the n x k matrix of posterior mixture assignment
// probabilities given current estimates of the model parameters for
// the special case when the residual variance is the same for all
// samples.
void compute_posterior_probs_iid (const mat& X, const vec& w, const cube& U,
				  const mat& V, mat& P) {
  
  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;

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
void compute_posterior_probs_notiid (const mat& X, const vec& w, 
				     const cube& U, const cube& V, 
				     mat& P) {
  
  // Get the number of rows (n) and columns (m) of the data matrix, and
  // the number of components in the mixture prior (k).
  unsigned int n = X.n_rows;
  unsigned int m = X.n_cols;
  unsigned int k = w.n_elem;

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
      P(i,j) = log(w(j)) + ldmvnorm(x,L);
    }
  }

  // Normalize the probabilities so that each row of P sums to 1.
  for (unsigned int i = 0; i < n; i++)
    P.row(i) = softmax(P.row(i));
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a mixture of
// multivariate normals, each with zero mean, covariance U[,,i] and
// weight w[i]. Return the posterior mean (mu1) and covariance (S1) of
// z. Note that input w1 must be the vector of *posterior* mixture
// weights (see compute_posterior_probs), and B1[,,i] should be the
// posterior covariance matrix for mixture component i. Input y is
// used to store intermediate calculations; it is a vector of the same
// size as x.
void compute_posterior_mvtnorm_mix (const vec& x, const vec& w1, 
				    const mat& V, const cube& B1,
				    vec& mu1, mat& S1) {
  unsigned int n = x.n_elem;
  unsigned int k = w1.n_elem;
  vec y(n);
  mu1.fill(0);
  S1.fill(0);
  for (unsigned int i = 0; i < k; i++) {
    compute_posterior_mean_mvtnorm(x,B1.slice(i),V,y);
    mu1 += w1(i) * y;
    S1  += w1(i) * (B1.slice(i) + y * trans(y));
  }
  S1 -= mu1 * trans(mu1);
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a multivariate normal
// distribution with mean zero and covariance U. Return in S1 the
// posterior covariance of z. (Note that the posterior covariance does
// not depend on x, so it is not one of the inputs.) These calculations 
// will only work if V is positive definite (invertible). Input I
// should be the identity matrix of the same dimension as U and V.
void compute_posterior_covariance_mvtnorm (const mat& U, const mat& V, 
					   mat& S1) {
  unsigned int n = U.n_rows;
  mat I(n,n,fill::eye);
  S1 = inv(U*inv(V) + I)*U;
}

// Suppose x is drawn from a multivariate normal distribution with
// mean z and covariance V, and z is drawn from a multivariate normal
// distribution with mean zero and covariance U. Given S1, the
// previously calculated posterior covariance of z, return the
// posterior mean (mu1) of z.
void compute_posterior_mean_mvtnorm (const vec& x, const mat& S1,
				     const mat& V, vec& mu1) {
  mu1 = S1 * solve(V,x);
}

