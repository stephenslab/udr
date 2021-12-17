#ifndef INCLUDE_MISC
#define INCLUDE_MISC

// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void         scale_rows    (arma::mat& A, const arma::vec& b);
double       logsumexp     (const arma::vec& x);
void         safenormalize (arma::vec& x);
arma::mat    crossprod     (const arma::mat& X);
double       ldmvnorm      (const arma::vec& x, const arma::mat& L);
double       ldmvnormmix   (const arma::vec& x, const arma::vec& w, 
			    const arma::cube& U, const arma::mat& V);
void         shrink_cov    (const arma::mat& T, arma::mat& U, 
			    double minval, unsigned int r);

#endif
