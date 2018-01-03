#ifndef LH_LOG_FUNCTIONS_H
#define LH_LOG_FUNCTIONS_H

#include <vector>
#include "row_set.hpp"

using namespace std;

double logdiff(double a, double b);

double logsum(double a, double b);

double log_big_sum(const vector<double>& R);

double log_big_sum(rowSet::const_iterator begin, rowSet::const_iterator end,
                   const vector<double>& R);

#endif