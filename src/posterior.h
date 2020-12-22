#ifndef INCLUDE_POSTERIOR
#define INCLUDE_POSTERIOR

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void compute_posterior_mvtnorm_mix (const arma::vec& x, const arma::vec& w1,
				    const arma::mat& V, const arma::cube& B1,
				    arma::vec& mu1, arma::mat& S1,
				    arma::vec& y);

void compute_posterior_covariance_mvtnorm (const arma::mat& U,
					   const arma::mat& V,
					   const arma::mat& I,
					   arma::mat& S1);

void compute_posterior_mean_mvtnorm (const arma::vec& x, const arma::mat& S1, 
				     const arma::mat& V, arma::vec& mu1);

#endif
