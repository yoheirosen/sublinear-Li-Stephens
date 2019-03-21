#include <cmath>
#include "math.hpp"
#include "penalty_set.hpp"

penaltySet::~penaltySet() {
  
}

penaltySet::penaltySet(double rho_in, double mu, int k) :
          k(k), mu(mu), rho(rho - log(k - 1)) {
  // basics
  rho = rho_in - log(k - 1);
  log_k = log(k);
  log_k_1 = log(k - 1);
  
  mu_c = log1p(-4*exp(mu));
  rho_c = log1p(-k*exp(rho));
  
  one_minus_mu = log1p(-4*exp(mu));
  one_minus_2mu = log1p(-5*exp(mu));
  row_coefficient = log1p(-k*exp(rho));
  rho_over_R_coeff = rho - row_coefficient;
  one_minus_mu_times_R_coeff = one_minus_mu + row_coefficient;
  mu_times_R_coeff = mu + row_coefficient;
}

double penaltySet::pow_rho_c(size_t l) const {
  return l * rho_c;
}

double penaltySet::pow_mu_c(size_t l) const {
  return l * mu_c;
}

double penaltySet::pow_mu(size_t l) const {
  return l * mu;
}

double penaltySet::span_polynomial(size_t l) const {
  return logmath::logdiff(0, pow_rho_c(l)) - log_k;
}

DPUpdateMap penaltySet::match_map(double last_sum) const {
  return DPUpdateMap(one_minus_mu_times_R_coeff, rho_over_R_coeff + last_sum);
}

DPUpdateMap penaltySet::non_match_map(double last_sum) const {
  return DPUpdateMap(mu_times_R_coeff, rho_over_R_coeff + last_sum);
}

DPUpdateMap penaltySet::current_map(double last_sum, bool match_is_rare) const {
  if(match_is_rare) {
    return non_match_map(last_sum);
  } else {
    return match_map(last_sum);
  }
}

double penaltySet::minority_correction(bool match_is_rare) const {
  if(match_is_rare) {
    return one_minus_mu - mu;
  } else {
    return mu - one_minus_mu;
  }
}

void penaltySet::update_sum(double& sum, const vector<double>& summands, 
              bool match_is_rare) const {
  if(match_is_rare) {
    double correct_to_1_m_2mu = one_minus_2mu - one_minus_mu;
    sum += mu;
    sum = logmath::logsum(sum, correct_to_1_m_2mu + logmath::log_big_sum(summands));
  } else {
    double correct_to_1_m_2mu = one_minus_2mu - mu;
    sum += one_minus_mu;
    sum = logmath::logdiff(sum, correct_to_1_m_2mu + logmath::log_big_sum(summands));
  }
}

double penaltySet::composed_row_coefficient(size_t l) const {
  return row_coefficient * l;
}

double penaltySet::span_mutation_penalty(size_t l) const {
  return l * one_minus_mu;
}

double penaltySet::span_coefficient(size_t l) const {
  return log1p(-exp(composed_row_coefficient(l))) - log_k;
}