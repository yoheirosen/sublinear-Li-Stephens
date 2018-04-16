#ifndef PENALTY_SET_H
#define PENALTY_SET_H

#include "math.hpp"
#include "DP_map.hpp"

using namespace std;

// stores a shared set of penalty-derived coefficients for use in calculations
// according to our model
struct penaltySet{
  int k;
  double rho;
  double rho_c;
  double mu;
  double mu_c;
  double log_k;
  double log_k_1;
  
  double pow_mu_c(size_t l) const;
  double pow_rho_c(size_t l) const;
  double pow_mu(size_t l) const;
  double span_polynomial(size_t l) const;
  
  double one_minus_mu;
  double one_minus_2mu;
  double rho_over_R_coeff;
  double one_minus_mu_times_R_coeff;
  double mu_times_R_coeff;
  
  double row_coefficient;
  
  penaltySet(double logRho, double logMu, int k);
  ~penaltySet();
  
  double composed_row_coefficient(size_t l) const;
  double span_mutation_penalty(size_t l) const;
  double span_coefficient(size_t l) const;
  
  DPUpdateMap match_map(double last_sum) const;
  DPUpdateMap non_match_map(double last_sum) const;
  DPUpdateMap current_map(double last_sum, bool match_is_rare) const;
  double minority_correction(bool match_is_rare) const;
  void update_sum(double& sum, const vector<double>& summands, bool match_is_rare) const;
  void update_sum(double& sum, const vector<double>& summands, rowSet::const_iterator begin, rowSet::const_iterator end, bool match_is_rare) const;
  
  // double mu_val(alleleValue from, alleleValue to) const;
  // double mu_loss_val(alleleValue from) const;
  // double rho_val(size_t position) const;
  // double rho_loss_val(size_t position) const;
};

#endif