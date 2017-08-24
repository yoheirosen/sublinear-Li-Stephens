#ifndef LH_DP_STATE_MAP
#define LH_DP_STATE_MAP

#include "lh_math.hpp"

struct DPUpdateMap{
private:
  bool degenerate_constant = false;
public:
  double coefficient;
  double constant;

  DPUpdateMap();
  DPUpdateMap(double coefficient);
  DPUpdateMap(double coefficient, double constant);
  DPUpdateMap(const DPUpdateMap& other);

  bool is_identity() const;
  bool is_degenerate() const;

  double of(double x) const;
  DPUpdateMap of(const DPUpdateMap& inner) const;

  // Composes the maps f1: x |-> A1(x + B1) and f2: x |-> A2(x + B2) to form a map
  // f': x |-> A2A1(x + B1 + B2/A1). Performs this in log-space
  DPUpdateMap compose(const DPUpdateMap& inner) const;
  void compose_in_place(const DPUpdateMap& inner);
  
  DPUpdateMap scale(double C) const;
  void scale_in_place(double C);
  
  bool operator==(const DPUpdateMap &other) const;
  bool operator!=(const DPUpdateMap &other) const;
};

#endif