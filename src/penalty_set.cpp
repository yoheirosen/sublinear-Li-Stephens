#include <cmath>
#include "penalty_set.hpp"

penaltySet::~penaltySet() {
  
}

penaltySet::penaltySet(double rho_in, double mu, int H) : H(H), 
          mu(mu) {
  rho = rho_in - log(H - 1);
  log_H = log(H);
  one_minus_mu = log1p(-4*exp(mu));
  one_minus_2mu = log1p(-5*exp(mu));
  R_coefficient = log1p(-H*exp(rho));
  rho_over_R_coeff = rho - R_coefficient;
  one_minus_mu_times_R_coeff = one_minus_mu + R_coefficient;
  mu_times_R_coeff = mu + R_coefficient;
}

DPUpdateMap penaltySet::get_match_map(double last_sum) const {
  return DPUpdateMap(one_minus_mu_times_R_coeff, rho_over_R_coeff + last_sum);
}

DPUpdateMap penaltySet::get_non_match_map(double last_sum) const {
  return DPUpdateMap(mu_times_R_coeff, rho_over_R_coeff + last_sum);
}

DPUpdateMap penaltySet::get_current_map(double last_sum, bool match_is_rare) const {
  if(match_is_rare) {
    return get_non_match_map(last_sum);
  } else {
    return get_match_map(last_sum);
  }
}

double penaltySet::get_minority_map_correction(bool match_is_rare) const {
  if(match_is_rare) {
    return one_minus_mu - mu;
  } else {
    return mu - one_minus_mu;
  }
}

void penaltySet::update_S(double& S, const vector<double>& summands, 
              bool match_is_rare) const {
  if(match_is_rare) {
    double correct_to_1_m_2mu = one_minus_2mu - one_minus_mu;
    S += mu;
    S = logsum(S, correct_to_1_m_2mu + log_big_sum(summands));
  } else {
    double correct_to_1_m_2mu = one_minus_2mu - mu;
    S += one_minus_mu;
    S = logdiff(S, correct_to_1_m_2mu + log_big_sum(summands));
  }
}

double penaltySet::composed_R_coefficient(size_t l) const {
  return R_coefficient * l;
}

double penaltySet::span_mutation_penalty(size_t l, size_t a) const {
  return (l - a) * one_minus_mu + a * mu;
}

double penaltySet::span_coefficient(size_t l) const {
  return log1p(-exp(composed_R_coefficient(l))) - log_H;
}