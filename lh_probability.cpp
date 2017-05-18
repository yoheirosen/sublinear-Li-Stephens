#include <cmath>
#include "lh_probability.hpp"

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

double log_big_sum(vector<double>& R) {
  if(R.size() == 0) {
    return nan("");
  } else if(R.size() == 1) {
    return R[0];
  } else {
    double max_summand = R[0];
    size_t max_index = 0;
    vector<double> summands;
    for(size_t i = 0; i < R.size(); i++){
      summands.push_back(R[i]);
      if(summands.back() > max_summand) {
        max_summand = summands.back();
        max_index = i;
      }
    }
    double sum = 0;
    for(size_t i = 0; i < summands.size(); i++) {
      if(i != max_index) {
        sum += exp(summands[i]-max_summand);
      }
    }
    return max_summand + log1p(sum);
  }
}

double log_weighted_big_sum(vector<double>& R, vector<size_t>& counts) {
  vector<double> new_R;
  for(size_t i = 0; i < R.size(); i++) {
    new_R.push_back(R[i] + log(counts[i]));
  }
  return log_big_sum(new_R);
}

void haplotypeMatrix::initialize_probability() {
  if(reference->has_span_before(0)) {
    int ell = reference->span_length_before(0);
    int augs = query->get_augmentations(-1);
    double fs = (ell - 1) * penalties->log_rho_complement;
    double ft = ell * penalties->log_ft_base;
    double mutation_coefficient = 
              (ell - augs) * (penalties->log_mu_complement) -
              augs * (penalties->log_mu);
    S[0] = mutation_coefficient + logsum(logsum(fs, ft), fs);
    double initial_value = mutation_coefficient +
              logsum(fs - penalties->log_H, penalties->log_rho + ft);
    initial_R = vector<double>(cohort->size(), initial_value);
    last_span_extended = -1;
  } else {
    initial_R = vector<double>(cohort->size(), 
              penalties->log_H - penalties->log_ft_base);
  }
}

void haplotypeMatrix::extend_probability_at_site(int j) {
  extend_probability_at_site(j, query->get_allele(j));
}

void haplotypeMatrix::extend_probability_at_site(int j, alleleValue a) {
  vector<size_t> matches = cohort->get_matches(j, a);
  // exp(C1) = (1-mu)(1-2rho)
  double C1 = penalties->log_R_match_coefficient;
  double C1m;
  // exp(C2) = (1-mu)rho * S
  double C2 = penalties->log_match_rho + last_S();
  double C2m;
  bool has_mismatches = (matches.size() != cohort->size());
  if(has_mismatches) {
    C1m = penalties->log_R_mismatch_coefficient;
    C2m = penalties->log_mismatch_rho + last_S();
  }
  for(size_t i = 0; matches.size(); i++) {
    R[j][i] = logsum(C1 + last_R(matches[i]),C2);
  }
  if(has_mismatches) {
    vector<int> mismatches;
    for(size_t i = 0; i < cohort->size(); i++) {
      if(R[j][i] == 0) {
        R[j][i] = logsum(C1m + last_R(matches[i]),C2m);
        mismatches.push_back(i);
      }
    }
    // Calculate S
    double log_M = log(matches.size());
    double log_Mm = log(mismatches.size());
    double C1_L = penalties->log_mu_complement +
            logsum(penalties->log_ft_base, penalties->log_rho + log_M);
    double C1_R = penalties->log_rho + penalties->log_mu + log_Mm;
    double C1 = last_S() + logdiff(C1_L,C1_R);
    double C2a = penalties->log_2mu_complement + penalties->log_ft_base;
    double C2b;
    vector<double> summands;
    // never sum more than 0.5* cohort-> size summands
    if(2 * mismatches.size() < cohort->size()) {
      for(size_t i = 0; i < mismatches.size(); i++) {
        summands.push_back(last_R(mismatches[i]));
      }
      C2b = log_big_sum(summands);
    } else {
      for(size_t i = 0; i < matches.size(); i++) {
        summands.push_back(last_R(matches[i]));
      }
      C2b = logdiff(last_S(), log_big_sum(summands));
    }
    S[j] = logdiff(C1, C2a + C2b);
  } else {
    // Calculate S, knowing that there are no mismatches
    S[j] = penalties->log_mu_complement + last_S() + penalties->log_fs_base;
  }
  last_extended = j;
  return;
}

void haplotypeMatrix::extend_probability_at_span_after(int j,
            int augmenation_count) {
  int ell = reference->span_length_after(j);
  double fs = (ell - 1) * penalties->log_fs_base;
  double ft = ell * penalties->log_ft_base;
  double C1 = (ell - augmenation_count) * penalties->log_mu_complement +
              augmenation_count * penalties->log_mu;
  double C2a = fs + penalties->log_rho + last_S();
  double C2b = penalties->log_ft_base + last_S();
  double C2 = logsum(C2a, C2b);  
  double diff = logdiff(penalties->log_ft_base + fs, ft);
  S[j] = C1 + logsum(C2, diff);
  double D1 = C1 + ft;
  double D2 = logsum(C2a, diff - penalties->log_H) - ft;
  for(size_t i = 0; i < cohort->size(); i++) {
    R[j][i] = D1 + logsum(R[j][i], D2);
  }
  last_span_extended = j;
}

bool haplotypeMatrix::last_extended_is_span() {
  return (last_extended == last_span_extended);
}

double haplotypeMatrix::last_S() {
  if(last_extended != -1) {
    return S[last_extended];
  } else {
    return S[0];
  }
}

double haplotypeMatrix::last_R(int i) {
  if(last_extended != -1) {
    return R[last_extended][i];
  } else {
    return initial_R[i];
  }
}

size_t haplotypeMatrix::size() {
  return reference->number_of_sites();
}

penaltySet::penaltySet(double log_rho, double log_mu, int H) : H(H), 
          log_rho(log_rho), log_mu(log_mu) {
  log_H = log(H);
  log_rho_complement = log1p(exp(log_rho));
  log_mu_complement = log1p(exp(log_mu));
  log_2mu_complement = log1p(2*exp(log_mu));
  log_ft_base = log1p(2*exp(log_rho));
  log_fs_base = logsum(log_ft_base, log_rho + log_H);
  log_R_match_coefficient = log_mu_complement + log_ft_base;
  log_R_mismatch_coefficient = log_mu + log_ft_base;
  log_match_rho = log_mu_complement + log_rho;
  log_mismatch_rho = log_mu + log_rho;
}

double haplotypeMatrix::calculate_probabilities() {
  initialize_probability();
  for(size_t i = 0; i < size(); i++) {
    extend_probability_at_site(i, query->get_allele(i));
    if(reference->has_span_after(i)) {
      extend_probability_at_span_after(i, query->get_augmentations(i));
    }
  }
  return S.back();
}

haplotypeMatrix::haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
          haplotypeCohort* cohort, inputHaplotype* query) :
          reference(ref), cohort(cohort), penalties(pen), query(query) {
  S = vector<double>(ref->number_of_sites(), 0);
  R = vector<vector<double> >(ref->number_of_sites(),
            vector<double>(cohort->size(),0));
}

// TODO: implement these 
// double haplotypeMatrix::minContinuing(int j) {
//   double newC = cohort->number_matching(j, query->get_allele(j));
//   double oldC = cohort->number_matching(j-1, query->get_allele(j-1));
//   return min(newC, oldC);
// }
// 
// double haplotypeMatrix::minMutating(int j) {
//   double newM = cohort->number_not_matching(j, query->get_allele(j));
//   double oldM = cohort->number_not_matching(j-1, query->get_allele(j-1));
//   return min(newM, oldM);
// }
// 
// double haplotypeMatrix::maxSwitching(int j) {
//   double newC = cohort->number_matching(j, query->get_allele(j));
//   double oldC = cohort->number_matching(j-1, query->get_allele(j-1));
//   return 0;
// }
// 
// double haplotypeMatrix::estimate_probability_at_site(int j, alleleValue a) {
//   return 0;
// }