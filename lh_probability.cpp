#include <iostream>
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

haplotypeMatrix::~haplotypeMatrix() {
  
}

void haplotypeMatrix::initialize_probability() {
  if(reference->has_span_before(0)) {
    // This span is initial: walks on the haplotype space begin with probability
    // |H|^-1 on a given haplotype rather than having their probabilities at the
    // left side of the span determined by history
    size_t length = reference->span_length_before(0);
    size_t num_augs = query->get_augmentations(-1);
    double l_fsl = (length - 1) * penalties->log_fs_base;
    double mutation_coefficient = 
              (length - num_augs) * (penalties->log_mu_complement) +
              num_augs * (penalties->log_mu);
    // Since we cannot distinguish haplotype prefixes over this initial span, we
    // know that all R-values at the right side of this span must be identical
    double initial_value = mutation_coefficient + l_fsl - penalties->log_H;
    initial_R = initial_value;
    S[0] = mutation_coefficient + l_fsl;
    last_span_extended = -1;
  }
}

void haplotypeMatrix::extend_probability_at_site(size_t j) {
  extend_probability_at_site(j, query->get_allele(j));
}

void haplotypeMatrix::extend_probability_at_site(size_t j, alleleValue a) {
  vector<size_t> matches = cohort->get_matches(j, a);
  bool has_mismatches = (matches.size() != cohort->size());
  double lm = penalties->log_mu;
  double lm_c = penalties->log_mu_complement;
  
  // Note that if there is an initial span, initialize_probability() sets the
  // counter last_span_extended to -1
  if(j == 0 && last_span_extended == -2) {
    // We are at an initial site; the probability of a (1-position-long) walk
    // begining and ending on a given haplotype is |H|^-1; the probability of it
    // emitting the query haplotype's allele at this site is simply governed by
    // whether a mutation is necessary
    
    // There are only two possible R-values at this site
    double initial_value = -penalties->log_H + lm_c;
    double m_initial_value = -penalties->log_H + lm;
    
    // Set the ones which match the query haplotype
    for(int i = 0; i < matches.size(); i++) {
      R[0][matches[i]] = initial_value;
    }
    // And the ones which don't
    for(int i = 0; i < cohort->size(); i++) {
      if(R[0][i] == 0) {
        R[0][i] = m_initial_value;
      }
    }
    
    // Get the sum of the R-values
    size_t num_mismatches = cohort->size() - matches.size();
    S[0] = -penalties->log_H + 
                logsum(log(matches.size()) + lm_c, log(num_mismatches) + lm);
    last_extended = 0;
    return;
  } else {
    // We are at a non-initial site, defined to include a site with index 0
    // which follows an initial span. 
    
    // This is the coefficient of the variable R[j-1][i] term
    double lft1 = penalties->log_ft_base;
    // This is the non-variable rho*S[j-1] term
    double lpS = penalties->log_rho + last_S();
    for(size_t i = 0; i < matches.size(); i++) {
      R[j][matches[i]] = lm_c + logsum(lft1 + last_R(matches[i]), lpS);
    }
    if(!has_mismatches) {
      // We updated all R-values already. Without any mismatches, S[j] takes a
      // very simple form (site j is just a span which got flagged as a site)
      S[j] = lm_c + last_S() + penalties->log_fs_base;
    } else {
      if(matches.size() == 0) {
        // We have updated no R-values. S[j] takes on a simple form in this case
        for(size_t i = 0; i < cohort->size(); i++) {
          R[j][i] = lm + logsum(lft1 + last_R(i), lpS);
        }
        S[j] = lm + last_S() + penalties->log_fs_base;
      } else {
        // This is a "real site." with variation in the populatino cohort. We 
        // need to find all the mismatch sites and update their R-values
        vector<size_t> mismatches;
        for(size_t i = 0; i < cohort->size(); i++) {
          if(R[j][i] == 0) {
            R[j][i] = lm + logsum(lft1 + last_R(i), lpS);
            mismatches.push_back(i);
          }
        }
        // We use an arithmetic trick to avoid ever summing more than 0.5*|H| 
        // terms. If there is a dominant allele in the population cohort at this
        // site we should be able to sum a very small number of terms!
        if(2 * mismatches.size() < cohort->size()) {
          // There are few mismatches; we want to sum them
          double log_num_mismatch = log(mismatches.size());
          double match_invariant = last_S() + penalties->log_fs_base +
                      penalties->log_mu_complement;
          double mismatch_invariant = log_num_mismatch + penalties->log_rho +
                      last_S();
          vector<double> summands;
          for(size_t i = 0; i < mismatches.size(); i++) {
            summands.push_back(last_R(mismatches[i]));
          }
          double mismatch_variant = lft1 + log_big_sum(summands);
          double corr_factor = logsum(mismatch_invariant, mismatch_variant) +
                      penalties->log_2mu_complement;
          S[j] = logdiff(match_invariant, corr_factor);
        } else {
          // There are few matches; we want to sum them
          double log_num_match = log(matches.size());
          double mismatch_invariant = last_S() + penalties->log_fs_base +
                      penalties->log_mu;
          double match_invariant = log_num_match + penalties->log_rho +
                      last_S();
          vector<double> summands;
          for(size_t i = 0; i < matches.size(); i++) {
            summands.push_back(last_R(matches[i]));
          }
          double match_variant = lft1 + log_big_sum(summands);
          double corr_factor = logsum(match_invariant, match_variant) +
                      penalties->log_2mu_complement;
          S[j] = logsum(mismatch_invariant, corr_factor);
        }
      }
    } 
    last_extended = j;
    return;
  }
}

void haplotypeMatrix::extend_probability_at_span_after(size_t j,
            int augmenation_count) {
  size_t length = reference->span_length_after(j);
  double l_fsl = length * penalties->log_fs_base;
  double l_ftl = length * penalties->log_ft_base;
  double mut_pen = (length - augmenation_count) * penalties->log_mu_complement +
              augmenation_count * penalties->log_mu;
  double l_sum_H = last_S() - penalties->log_H;
  double R_invariant = l_sum_H + logdiff(l_fsl, l_ftl);
  for(size_t i = 0; i < cohort->size(); i++) {
    R[j][i] = mut_pen + logsum(l_ftl + R[j][i], R_invariant);
  }
  S[j] = mut_pen + last_S() + l_fsl;
  last_span_extended = j;
}

bool haplotypeMatrix::last_extended_is_span() {
  return (last_extended == last_span_extended);
}

double haplotypeMatrix::last_S() {
  if(last_extended != -1) {
    return S[last_extended];
  } else {
    return initial_R + penalties->log_H;
  }
}

double haplotypeMatrix::last_R(size_t index_among_haplotypes) {
  if(last_extended != -1) {
    return R[last_extended][index_among_haplotypes];
  } else {
    return initial_R;
  }
}

size_t haplotypeMatrix::size() {
  return reference->number_of_sites();
}

penaltySet::~penaltySet() {
  
}

penaltySet::penaltySet(double log_rho, double log_mu, int H) : H(H), 
          log_rho(log_rho), log_mu(log_mu) {
  log_H = log(H);
  log_rho_complement = log1p(-exp(log_rho));
  log_mu_complement = log1p(-exp(log_mu));
  log_2mu_complement = log1p(-2*exp(log_mu));
  log_ft_base = log1p(-2*exp(log_rho));
  log_fs_base = logsum(log_ft_base, log_rho + log_H);
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
// double haplotypeMatrix::estimate_probability_at_site(size_t j, alleleValue a) {
//   return 0;
// }

// vector<size_t> get_matches(size_t index, alleleValue a) {
//   return cohort->get_matches(query->get_global_index(i), a);
// }

// size_t number_matching(size_t index, alleleValue a) {
//   return cohort->number_matching(query->get_global_index(i), a);
// }

// size_t number_not_matching(size_t index, alleleValue a) {
//   return cohort->number_not_matching(query->get_global_index(i), a);
// }