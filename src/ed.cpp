#include "misc.h"

using namespace Rcpp;
using namespace arma;

// FUNCTION DECLARATIONS
// ---------------------
void ed_iid (const mat& X, mat& U, const vec& p);

// FUNCTION DEFINITIONS
// --------------------
// Update the prior covariance matrix (U) in the model x ~ N(0,U + I)
// using the update formula derived in Bovy et al (2011). Input p is a
// vector of weights associated with the rows of X.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ed_iid_rcpp (const arma::mat& X, const arma::mat& U, 
		       const arma::vec& p) {
  mat Unew = U;
  ed_iid(X,Unew,p);
  return Unew;
}

// Update the prior covariance matrix (U) in the model x ~ N(0,U + I)
// using the update formula derived in Bovy et al (2011). Input p is a
// vector of weights associated with the rows of X.
void ed_iid (const mat& X, mat& U, const vec& p) {
  unsigned int m = X.n_cols;
  vec p1 = p;
  mat X1 = X;
  safenormalize(p1);
  scale_rows(X1,sqrt(p1));
  mat I(m,m,fill::eye);
  mat T = U + I;
  mat B = solve(T,U);
  X1 *= B;
  U += crossprod(X1) - U*B;
}
