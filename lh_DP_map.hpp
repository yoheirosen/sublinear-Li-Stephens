#ifndef LH_DP_STATE_MAP
#define LH_DP_STATE_MAP

#include "lh_math.hpp"

struct DPUpdateMap{
  double coefficient;
  double constant;

  DPUpdateMap();
  DPUpdateMap(double coefficient, double constant);

  double evaluate_at(double x);

  // Composes the maps f1: x |-> A1(x + B1) and f2: x |-> A2(x + B2) to form a map
  // f': x |-> A2A1(x + B1 + B2/A1). Performs this in log-space
  DPUpdateMap compose(DPUpdateMap inner);
};

#endif