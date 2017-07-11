#include "lh_math.hpp"
#include "lh_DP_map.hpp"

DPUpdateMap::DPUpdateMap() {}

DPUpdateMap::DPUpdateMap(double coefficient, double constant) : 
          coefficient(coefficient), constant(constant) {  
}

double DPUpdateMap::evaluate_at(double x) {
  return coefficient + logsum(x, constant);
}

// Composes the maps f1: x |-> A1(x + B1) and f2: x |-> A2(x + B2) to form a map
// f': x |-> A2A1(x + B1 + B2/A1). Performs this in log-space
DPUpdateMap DPUpdateMap::compose(DPUpdateMap outer, DPUpdateMap inner) {
  DPUpdateMap to_return;
  to_return.coefficient = outer.coefficient + inner.coefficient; 
  to_return.constant = logsum(inner.constant, 
                              outer.constant - inner.coefficient);
  return to_return;
}