#include <cmath>
#include "lh_probability.hpp"

haplotypeMatrix::~haplotypeMatrix() {
  
}

void haplotypeMatrix::initialize_probability() {
  if(query->has_left_tail()) {
    // This span is initial: walks on the haplotype space begin with probability
    // |H|^-1 on a given haplotype rather than having their probabilities at the
    // left side of the span determined by history
    
    size_t length = query->get_left_tail();
    size_t num_augs = query->get_augmentations(-1);
    double lfsl = (length - 1) * penalties->log_fs_base;
    double mutation_coefficient = 
              (length - num_augs) * (penalties->log_mu_complement) +
              num_augs * (penalties->log_mu);
    // Since we cannot distinguish haplotype prefixes over this initial span, we
    // know that all R-values at the right side of this span must be identical
    double initial_value = mutation_coefficient + lfsl - penalties->log_H;
    initial_R = initial_value;
    S[0] = mutation_coefficient + lfsl;
    last_span_extended = -1;
  }
}

void haplotypeMatrix::extend_probability_at_site(size_t j) {
  extend_probability_at_site(j, query->get_allele(j));
}

double calculate_R(double oldR, DPUpdateMap map) {
  return map.evaluate_at(oldR);
}

double calculate_R(double oldR, double coefficient, double constant) {
  return calculate_R(oldR, DPUpdateMap(coefficient, constant));
}

void haplotypeMatrix::extend_probability_at_site(size_t j, alleleValue a) {
  double lm = penalties->log_mu;
  double lm_c = penalties->log_mu_complement;

  if(j == 0 && last_span_extended == -2) {
    vector<size_t> matches = get_matches(j);
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
    size_t num_mismatches = cohort->number_not_matching(j, a);
    S[0] = -penalties->log_H + 
                logsum(log(matches.size()) + lm_c, log(num_mismatches) + lm);
    last_extended = 0;
    return;
  } else {
    // This is the coefficient of the variable R[j-1][i] term
    double lft1 = penalties->log_ft_base;
    // This is the non-variable rho*S[j-1] term
    double lpS = penalties->log_rho + last_S();    
    if(cohort->match_is_rare(j, a)) {
      map.add_map_for_site(lm + lft1, lpS - lft1);
      if(cohort->number_matching(j, a) == 0) {
        S[j] = lm + last_S() + penalties->log_fs_base;
      } else {
        vector<size_t> matches = cohort->get_matches(j, a);
        vector<double> summands;
        for(size_t i = 0; i < matches.size(); i++) {
          summands.push_back(last_R(matches[i]));
        }
        
        vector<size_t> slots = map.rows_to_slots(matches);
        map.update_maps(slots);
        // it is important to note that at this point, the coefficient of the maps
        // is off by a factor of lm_c - lm for all matching rows
        double flip_m = lm_c - lm;
        for(size_t i = 0; i < matches.size(); i++) {
          R[j][matches[i]] = flip_m +
                      calculate_R(last_R(matches[i]), map.get_map(matches[i]));        
          map.remove_row_from_slot(matches[i]);
        }
        map.add_identity_map();
        for(size_t i = 0; i < matches.size(); i++) {
          map.assign_row_to_newest_index(matches[i]);
        }
        
        double log_num_match = log(cohort->number_matching(j,a));
        double mismatch_invariant = last_S() + penalties->log_fs_base +
                    penalties->log_mu;
        double match_invariant = log_num_match + penalties->log_rho +
                    last_S();
        double match_variant = lft1 + log_big_sum(summands);
        double corr_factor = logsum(match_invariant, match_variant) +
                    penalties->log_2mu_complement;
        S[j] = logsum(mismatch_invariant, corr_factor);
      }
    } else {
      map.add_map_for_site(lm_c + lft1, lpS - lft1);
      if(cohort->number_not_matching(j, a) == 0) {
        S[j] = lm_c + last_S() + penalties->log_fs_base;
      } else {
        vector<size_t> non_matches = cohort->get_non_matches(j, a);
        vector<double> summands;
        for(size_t i = 0; i < non_matches.size(); i++) {
          summands.push_back(last_R(non_matches[i]));
        }
        
        vector<size_t> slots = map.rows_to_slots(non_matches);
        map.update_maps(slots);
        double flip_m = lm - lm_c;
        for(size_t i = 0; i < non_matches.size(); i++) {
          R[j][non_matches[i]] = flip_m +
                    calculate_R(last_R(non_matches[i]), 
                                map.get_map(non_matches[i]));        
          map.remove_row_from_slot(non_matches[i]);
        }
        map.add_identity_map();
        for(size_t i = 0; i < non_matches.size(); i++) {
          map.assign_row_to_newest_index(non_matches[i]);
        }
      
        double log_num_mismatch = log(cohort->number_not_matching(j,a));
        double match_invariant = last_S() + penalties->log_fs_base +
                    penalties->log_mu_complement;
        double mismatch_invariant = log_num_mismatch + penalties->log_rho +
                    last_S();
        double mismatch_variant = lft1 + log_big_sum(summands);
        double corr_factor = logsum(mismatch_invariant, mismatch_variant) +
                    penalties->log_2mu_complement;
        S[j] = logdiff(match_invariant, corr_factor);      
      }
    }
  }
  last_extended = j;
  return;
}

void haplotypeMatrix::extend_probability_at_span_after(size_t j,
            int augmenation_count) {
  size_t length = query->get_span_after(j);
  double lfsl = length * penalties->log_fs_base;
  double lftl = length * penalties->log_ft_base;
  double mut_pen = (length - augmenation_count) * penalties->log_mu_complement +
              augmenation_count * penalties->log_mu;
  double lsum_H = last_S() - penalties->log_H;
  double R_invariant = lsum_H + logdiff(lfsl, lftl);
  map.update_map_with_span(mut_pen + lftl, R_invariant - lftl);
  S[j] = mut_pen + last_S() + lfsl;
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
  return query->number_of_sites();
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

double haplotypeMatrix::calculate_probability() {
  initialize_probability();
  if(query->has_sites()) {
    for(size_t i = 0; i < size(); i++) {
      extend_probability_at_site(i, query->get_allele(i));
      if(query->has_span_after(i)) {
        extend_probability_at_span_after(i, query->get_augmentations(i));
      }
    }
  }
  return S.back();
}

haplotypeMatrix::haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
          haplotypeCohort* cohort, inputHaplotype* query) :
          reference(ref), cohort(cohort), penalties(pen), query(query),
          map(delayMap(cohort->size(), 0)) {
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

vector<size_t> haplotypeMatrix::get_matches(size_t i) {
  return cohort->get_matches(query->get_rel_index(i), query->get_allele(i));
}

size_t haplotypeMatrix::number_matching(size_t i) {
  return cohort->number_matching(query->get_rel_index(i), query->get_allele(i));
}

size_t haplotypeMatrix::number_not_matching(size_t i) {
  return cohort->number_not_matching(query->get_rel_index(i),
            query->get_allele(i));
}