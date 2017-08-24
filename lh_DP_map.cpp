#include "lh_math.hpp"
#include "lh_DP_map.hpp"

DPUpdateMap::DPUpdateMap() {}

DPUpdateMap::DPUpdateMap(double coefficient) : coefficient(coefficient) {
  degenerate_constant = true;
  constant == 0;
}

DPUpdateMap::DPUpdateMap(double coefficient, double constant) : 
          coefficient(coefficient), constant(constant) {
}

DPUpdateMap::DPUpdateMap(const DPUpdateMap& other) {
  degenerate_constant = other.degenerate_constant;
  coefficient = other.coefficient;
  constant = other.constant;
}

double DPUpdateMap::of(double x) const {
  if(degenerate_constant) {
    return coefficient + x;
  } else {
    return coefficient + logsum(x, constant);
  }
}

DPUpdateMap DPUpdateMap::of(const DPUpdateMap& inner) const {
  DPUpdateMap to_return = *this;
  to_return.compose_in_place(inner);
  return to_return;
}

DPUpdateMap DPUpdateMap::compose(const DPUpdateMap& inner) const {
  return this->of(inner);
}

void DPUpdateMap::compose_in_place(const DPUpdateMap& inner) {
  if(degenerate_constant && inner.degenerate_constant) {
    coefficient = coefficient + inner.coefficient;
    return;
  } else if(degenerate_constant) {
    coefficient = coefficient + inner.coefficient;
    constant = inner.constant;  
    degenerate_constant = false;
    return;
  } else if(inner.degenerate_constant) {
    coefficient = coefficient + inner.coefficient;
    constant = constant - inner.coefficient;
    return;
  } else {
    coefficient = coefficient + inner.coefficient; 
    constant = logsum(inner.constant, constant - inner.coefficient);
    return;
  }
}

DPUpdateMap DPUpdateMap::scale(double C) const {
  DPUpdateMap to_return = *this;
  to_return.scale_in_place(C);
  return to_return;
}

void DPUpdateMap::scale_in_place(double C) {
  coefficient += C;
}

bool DPUpdateMap::is_identity() const {
  return coefficient == 0 && degenerate_constant;
}

bool DPUpdateMap::is_degenerate() const {
  return degenerate_constant;
}


bool DPUpdateMap::operator==(const DPUpdateMap &other) const {
  if(degenerate_constant && other.degenerate_constant) {
    return coefficient == other.coefficient;
  } if(degenerate_constant != other.degenerate_constant) {
    return false;
  } else {
    return constant == other.constant && coefficient == other.coefficient;
  }
}

bool DPUpdateMap::operator!=(const DPUpdateMap &other) const {
  return !(*this == other);
}