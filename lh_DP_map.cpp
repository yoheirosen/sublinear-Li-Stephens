#include "lh_math.hpp"
#include "lh_DP_map.hpp"

DPUpdateMap::DPUpdateMap() {}

DPUpdateMap::DPUpdateMap(double coefficient, double constant) : 
          coefficient(coefficient), constant(constant) {  
}

double DPUpdateMap::evaluate_at(double x) {
  if(constant < 0 && coefficient == 0) {
    return x;
  } else {
    return coefficient + logsum(x, constant);
  }
}

DPUpdateMap DPUpdateMap::compose(DPUpdateMap inner) {
  if(constant < 0 && coefficient == 0) {
    return inner;
  } else if(inner.constant < 0 && inner.coefficient == 0){
    return *this;
  } else {
    DPUpdateMap to_return;
    to_return.coefficient = coefficient + inner.coefficient; 
    to_return.constant = logsum(constant, 
                                constant - inner.coefficient);
    return to_return;
  }
}