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

DPUpdateMap DPUpdateMap::compose(DPUpdateMap& inner) {
  DPUpdateMap to_return = *this;
  to_return.compose_in_place(inner);
  return to_return;
}

void DPUpdateMap::compose_in_place(DPUpdateMap& inner) {
  if(constant < 0 && coefficient == 0) {
    *this = inner;
    return;
  } else if(inner.constant < 0 && inner.coefficient == 0){
    return;
  } else if(constant < 0) {
    double C = coefficient;
    *this = inner;
    this->scale_in_place(C);
    return;
  } else if(inner.constant < 0) {
    this->scale_in_place(coefficient);
    return;
  } else {
    this->coefficient = this->coefficient + inner.coefficient; 
    this->constant = logsum(inner.constant, 
                                this->constant - inner.coefficient);
    return;
  }
}

DPUpdateMap DPUpdateMap::scale(double C) {
  DPUpdateMap to_return = *this;
  to_return.scale_in_place(C);
  return to_return;
}

void DPUpdateMap::scale_in_place(double C) {
  coefficient += C;
}