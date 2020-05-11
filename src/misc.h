#ifndef INCLUDE_MISC
#define INCLUDE_MISC

#include <RcppArmadillo.h>

// FUNCTION DECLARATIONS
// ---------------------
void         scale_rows (arma::mat& A, const arma::vec& b);
arma::rowvec softmax    (const arma::rowvec & x);
double       ldmvnorm   (const arma::vec& x, const arma::mat& L);

#endif
