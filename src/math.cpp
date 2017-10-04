#include <cmath>
#include <vector>
#include "math.hpp"

using namespace std;

double logdiff(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(-exp(b - a));
}

double logsum(double a, double b) {
  if(b > a) {
    double c = a;
    a = b;
    b = c;
  }
  return a + log1p(exp(b - a));
}

double log_big_sum(const vector<double>& R) {
  if(R.size() == 0) {
    return nan("");
  } else if(R.size() == 1) {
    return R[0];
  } else {
    double max_summand = R[0];
    size_t max_index = 0;
    for(size_t i = 0; i < R.size(); i++){
      if(R[i] > max_summand) {
        max_summand = R[i];
        max_index = i;
      }
    }
    double sum = 0;
    for(size_t i = 0; i < R.size(); i++) {
      if(i != max_index) {
        sum += exp(R[i] - max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}

double log_weighted_big_sum(const vector<double>& R, const vector<size_t>& counts) {
  vector<double> new_R;
  for(size_t i = 0; i < R.size(); i++) {
    new_R.push_back(R[i] + log(counts[i]));
  }
  return log_big_sum(new_R);
}
