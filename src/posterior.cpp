#include "misc.h"
#include "posterior.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
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

