#include "misc.h"

using namespace arma;

// FUNCTION DEFINITIONS
// --------------------
// Return the "softmax" of vector x, y(i) = exp(x(i))/sum(exp(x)), in
// a way that guards against numerical underflow or overflow. The
// return value is a vector with entries that sum to 1.
vec softmax (const vec & x) {
  vec y = exp(x - max(x));
  y /= sum(y);
  return y;
}

