#include <cmath>
#include "penalty_set.hpp"

penaltySet::~penaltySet() {
  
}

penaltySet::penaltySet(double rho, double mu, int H) : H(H), 
          rho(rho), mu(mu) {
  log_H = log(H);
  one_minus_rho = log1p(-exp(rho));
  one_minus_mu = log1p(-exp(mu));
  one_minus_2mu = log1p(-2*exp(mu));
  alpha_value = log1p(-2*exp(rho));
  beta_value = logsum(alpha_value, rho + log_H);
}

DPUpdateMap penaltySet::get_match_map(double last_sum) const {
  return DPUpdateMap(one_minus_mu + alpha_value, rho + last_sum - alpha_value);
}

DPUpdateMap penaltySet::get_non_match_map(double last_sum) const {
  return DPUpdateMap(mu + alpha_value, rho + last_sum - alpha_value);
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
    S += beta_value + mu;
    S = logsum(S, correct_to_1_m_2mu + log_big_sum(summands));
  } else {
    double correct_to_1_m_2mu = one_minus_2mu - mu;
    S += beta_value + one_minus_mu;
    S = logdiff(S, correct_to_1_m_2mu + log_big_sum(summands));
  }
}