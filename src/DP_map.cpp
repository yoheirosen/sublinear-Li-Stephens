#include "math.hpp"
#include "DP_map.hpp"

DPUpdateMap::DPUpdateMap() {}

DPUpdateMap::DPUpdateMap(double coefficient) : coefficient(coefficient) {
  scalar = true;
  constant == 0;
}

DPUpdateMap::DPUpdateMap(double coefficient, double constant) : 
          coefficient(coefficient), constant(constant) {
}

DPUpdateMap::DPUpdateMap(const DPUpdateMap& other) {
  scalar = other.scalar;
  coefficient = other.coefficient;
  constant = other.constant;
}

double DPUpdateMap::of(double x) const {
  if(scalar) {
    return coefficient + x;
  } else {
    return coefficient + logmath::logsum(x, constant);
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
  if(scalar && inner.scalar) {
    coefficient = coefficient + inner.coefficient;
    return;
  } else if(scalar) {
    coefficient = coefficient + inner.coefficient;
    constant = inner.constant;  
    scalar = false;
    return;
  } else if(inner.scalar) {
    coefficient = coefficient + inner.coefficient;
    constant = constant - inner.coefficient;
    return;
  } else {
    coefficient = coefficient + inner.coefficient; 
    constant = logmath::logsum(inner.constant, constant - inner.coefficient);
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
  return coefficient == 0 && scalar;
}

bool DPUpdateMap::is_degenerate() const {
  return scalar;
}

bool DPUpdateMap::operator==(const DPUpdateMap &other) const {
  if(scalar && other.scalar) {
    return coefficient == other.coefficient;
  } if(scalar != other.scalar) {
    return false;
  } else {
    return constant == other.constant && coefficient == other.coefficient;
  }
}

bool DPUpdateMap::operator!=(const DPUpdateMap &other) const {
  return !(*this == other);
}

const DPUpdateMap DPUpdateMap::IDENTITY = DPUpdateMap(0);

// 
// DPUpdateMap& operator+=(const DPUpdateMap& other) {
//   
// }
// 
// DPUpdateMap operator+(const DPUpdateMap& other) const {
//   DPUpdateMap to_return = *this;
//   to_return += other;
//   return to_return;
// }
// 
// DPUpdateMap& operator-=(const DPUpdateMap& other);
// DPUpdateMap operator-(const DPUpdateMap& other) const;