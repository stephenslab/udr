#include "misc.h"

using namespace arma;

// INLINE FUNCTION DEFINITIONS
// ---------------------------
// Return a or b, which ever is smaller.
inline double min (double a, double b) {
  double y;
  if (a < b)
    y = a;
  else
    y = b;
  return y;
}

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
// Scale each row A[i,] by b[i].
void scale_rows (mat& A, const vec& b) {
  vec c = b;
  A.each_col() %= c;
}

// Compute log(sum(exp(x))) in a numerically stable way.
double logsumexp (const vec& x) {
  double a = max(x);
  return log(sum(exp(x - a))) + a;
}

// Replace x with x/sum(x), but take care of the special case when all
// the entries are zero, in which case return the vector of all 1/n,
// where n is the length of x.
void safenormalize (vec& x) {
  unsigned int n = x.n_elem;
  if (sum(x) <= 0)
    x.fill(1/n);
  else
    x = x/sum(x);
}

// Return the cross-product of matrix X, i.e., X'*X.
mat crossprod (const mat& X) {
  return trans(X) * X;
}

// Compute the log-probability of x, where x is multivariate normal
// with mean zero and covariance matrix S. Input argument should be L
// be the Cholesky factor of S; L = chol(S,"lower").
double ldmvnorm (const vec& x, const mat& L) {
  double d = norm(solve(L,x),2);
  return -d*d/2 - sum(log(sqrt(2*M_PI)*L.diag()));
}

// Compute mixture of normals log-density in a numerically stable
// way. Specifically, this function should return the same result as
//
//  k <- length(w)
//  log(sum(sapply(1:k,function (i) w[i] * dmvnorm(x,sigma = U[,,i] + V))))
//
// except that the computation is done in a more numerically stable
// way.
double ldmvnormmix (const vec& x, const vec& w, const cube& U, const mat& V) {
  unsigned int n = w.n_elem;
  unsigned int m = x.n_elem;
  vec y(n);
  mat T(m,m);
  mat L(m,m);
  for (unsigned int i = 0; i < n; i++) {
    T    = U.slice(i) + V;
    L    = chol(T,"lower");
    y(i) = log(w(i)) + ldmvnorm(x,L);
  }
  return logsumexp(y);
}

// Find the n x n matrix U + I that best approximates T satisfying the
// constraint that U is positive definite. This is achieved by setting
// any eigenvalues of T less than 1 to 1 + minval, or, equivalently,
// setting any eigenvalues of U less than zero to be minval. The output
// is a positive definite matrix, U. When r < n, U is additionally
// constrained so that at most r of its eigenvalues are positive.
void shrink_cov (const mat& T, mat& U, double minval, unsigned int r) {
  unsigned int n = T.n_rows;
  r = min(r,n);
  mat Y(n,n);
  vec d(n);
  eig_sym(d,Y,T);
  for (unsigned int i = 0; i < n; i++)
    d(i) = max(d(i) - 1,minval);
  if (r < n)
    for (unsigned int i = n - r; i-- > 0; ) 
      d(i) = minval;

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
