#ifndef LH_LOG_FUNCTIONS_H
#define LH_LOG_FUNCTIONS_H

#include <vector>
#include <cmath>

namespace logmath{

inline double logdiff(double a, double b) {
  if(b > a) { std::swap(a, b); }
  return a + log1p(-exp(b - a));
}

inline double logsum(double a, double b) {
  if(b > a) { std::swap(a, b); }
  return a + log1p(exp(b - a));
}

template<class IteratorType>
double log_big_sum(const IteratorType& begin, const IteratorType& end,
                   const std::vector<double>& R) {
  IteratorType it = begin;
  ++it;
  if(it == end) {
    return R[*begin];
  } else {
    double max_summand = R[*begin];
    IteratorType max_index = begin;
    for(it = begin; it != end; ++it){
      if(R[*it] > max_summand) {
        max_summand = R[*it];
        max_index = it;
      }
    }
    double sum = 0;
    for(it = begin; it != end; ++it) {
      if(it != max_index) {
        sum += exp(R[*it] - max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}

inline double log_big_sum(const std::vector<double>& R) {
  return log_big_sum(R.begin(), R.end(), R);
}

} // namespace log_math

#endif
