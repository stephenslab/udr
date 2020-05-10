#ifndef INCLUDE_LIKELIHOOD
#define INCLUDE_LIKELIHOOD

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
double ldmvnorm (const arma::vec& x, const arma::mat& L);

double loglik_mvebnm (const arma::mat& X, const arma::vec& w, 
		      const arma::cube& U, const arma::mat& S);

#endif
