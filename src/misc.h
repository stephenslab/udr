#ifndef INCLUDE_MISC
#define INCLUDE_MISC

// This is included to suppress the warnings from solve() when the
// system is singular or close to singular.
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void         scale_rows (arma::mat& A, const arma::vec& b);
arma::rowvec softmax    (const arma::rowvec & x);
arma::mat    crossprod  (const arma::mat& X);
double       ldmvnorm   (const arma::vec& x, const arma::mat& L);

#endif
