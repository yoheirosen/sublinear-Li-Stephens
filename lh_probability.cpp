#include <cmath>
#include "lh_probability.hpp"
#include <iostream>

haplotypeMatrix::~haplotypeMatrix() {
  
}

delayMap haplotypeMatrix::get_map() {
  return map;
}

void haplotypeMatrix::account_for_initial_span() {
  if(query->has_left_tail()) {
    // There is a uniform 1/|H| probability of starting on any given haplotype.
    // All emission probabilities are the same. So all R-values are the same
    
    size_t length = query->get_left_tail();
    size_t num_augs = query->get_augmentations(-1);
    double lfsl = (length - 1) * penalties->log_fs_base;
    double mutation_coefficient = 
              (length - num_augs) * (penalties->log_mu_complement) +
              num_augs * (penalties->log_mu);
    initial_R = mutation_coefficient + lfsl - penalties->log_H;
    for(size_t i = 0; i < R.size(); i++) {
      R[i] = initial_R;
    }
    S = mutation_coefficient + lfsl;
    last_span_extended = -1;
  }
}

void haplotypeMatrix::extend_first_site() {
  if(query->has_left_tail()) {
    extend_probability_at_site(0);
  } else {
    vector<size_t> matches = cohort->get_matches(0, query->get_allele(0));
    vector<size_t> non_matches = 
              cohort->get_non_matches(0, query->get_allele(0));
    // There are only two possible R-values at this site. There is a uniform
    // 1/|H| probability of starting on any given haplotype; the emission
    // probabilities account for differences in R-value
    double initial_value = -penalties->log_H + penalties->log_mu_complement;
    double m_initial_value = -penalties->log_H + penalties->log_mu;
    
    // Set the ones which match the query haplotype
    for(size_t i = 0; i < matches.size(); i++) {
      R[matches[i]] = initial_value;
    }
    // And the ones which don't
    for(size_t i = 0; i < non_matches.size(); i++) {
      R[non_matches[i]] = m_initial_value;
    }
    
    S = -penalties->log_H + 
                logsum(log(matches.size()) + penalties->log_mu_complement,
                       log(non_matches.size()) + penalties->log_mu);
    last_extended = 0;
  }
}


void haplotypeMatrix::initialize_probability() {
  account_for_initial_span();
  extend_first_site();
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

void haplotypeMatrix::take_snapshot() {
  map.hard_update_all();
  size_t j = last_extended;
  alleleValue a = query->get_allele(last_extended);
  bool homogenous = (cohort->number_matching(j, a) == 0 ||
            cohort->number_not_matching(j, a) == 0);
  if(homogenous || last_extended_is_span()) {
    for(int i = 0; i < cohort->size(); i++) {
      R[i] = calculate_R(R[i], map.get_map(i));
    }
  } else if(cohort->match_is_rare(j, a)) {
    vector<size_t> non_matches = cohort->get_non_matches(j, a);
    for(int i = 0; i < non_matches.size(); i++) {
      R[non_matches[i]] =
                calculate_R(R[non_matches[i]], map.get_map(non_matches[i]));
    }
  } else {
    vector<size_t> matches = cohort->get_matches(j, a);
    for(int i = 0; i < matches.size(); i++) {
      R[matches[i]] =
                calculate_R(R[matches[i]], map.get_map(matches[i]));
    }
  }
  map.hard_clear_all();
} 

void haplotypeMatrix::extend_probability_at_site(size_t j, alleleValue a) {
  double lm = penalties->log_mu;
  double lm_c = penalties->log_mu_complement;
  // This is the coefficient of the variable R[j-1][i] term
  double lft1 = penalties->log_ft_base;
  // This is the non-variable rho*S[j-1] term
  double lpS = penalties->log_rho + S;    

  if(cohort->match_is_rare(j, a)) {
    map.add_map_for_site(lm + lft1, lpS - lft1);
    if(cohort->number_matching(j, a) == 0) {
      // separate case to avoid log-summing "log 0"
      S = lm + S + penalties->log_fs_base;
    } else {
      vector<size_t> matches = cohort->get_matches(j, a);
            
      vector<size_t> slots = map.rows_to_slots(matches);
      map.update_maps(slots);
      // it is important to note that at this point, the coefficient of the maps
      // is off by a factor of lm_c - lm for all matching rows
      double flip_m = lm_c - lm;
      for(size_t i = 0; i < matches.size(); i++) {
        R[matches[i]] = flip_m +
                    calculate_R(R[matches[i]], map.get_map(matches[i]));        
        map.remove_row_from_slot(matches[i]);
      }
      map.add_identity_map();
      for(size_t i = 0; i < matches.size(); i++) {
        map.assign_row_to_newest_index(matches[i]);
      }
      
      vector<double> summands;
      for(size_t i = 0; i < matches.size(); i++) {
        summands.push_back(R[matches[i]]);
      }
      
      double l_num_match = log(cohort->number_matching(j,a));
      double l_2mc = penalties->log_2mu_complement;
      double l_2mc_flip = l_2mc - lm_c;
      double l_p = penalties->log_rho;
      
      S += penalties->log_fs_base + lm;
      S = logsum(S, l_2mc_flip + log_big_sum(summands));
    }
  } else {
    map.add_map_for_site(lm_c + lft1, lpS - lft1);
    if(cohort->number_not_matching(j, a) == 0) {
      // separate case to avoid log-summing "log 0"
      S = lm_c + S + penalties->log_fs_base;
    } else {
      vector<size_t> non_matches = cohort->get_non_matches(j, a);
      
      vector<size_t> slots = map.rows_to_slots(non_matches);
      map.update_maps(slots);
      double flip_m = lm - lm_c;
      for(size_t i = 0; i < non_matches.size(); i++) {
        R[non_matches[i]] = flip_m +
                  calculate_R(R[non_matches[i]], 
                              map.get_map(non_matches[i]));        
        map.remove_row_from_slot(non_matches[i]);
      }
      map.add_identity_map();
      for(size_t i = 0; i < non_matches.size(); i++) {
        map.assign_row_to_newest_index(non_matches[i]);
      }
    
      vector<double> summands;
      for(size_t i = 0; i < non_matches.size(); i++) {
        summands.push_back(R[non_matches[i]]);
      }
    
      double l_num_mismatch = log(cohort->number_not_matching(j,a));
      double l_2mc = penalties->log_2mu_complement;
      double l_2mc_flip = l_2mc - lm;
      double l_p = penalties->log_rho;
      
      S += penalties->log_fs_base + lm_c;
      S = logdiff(S, l_2mc_flip + log_big_sum(summands));     
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
  double RHS = S - penalties->log_H + logdiff(lfsl, lftl);
  map.update_map_with_span(mut_pen + lftl, RHS - lftl);
  S = mut_pen + S + lfsl;
  last_span_extended = j;
}

bool haplotypeMatrix::last_extended_is_span() {
  return (last_extended == last_span_extended);
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
    if(query->has_span_after(0)) {
      extend_probability_at_span_after(0, query->get_augmentations(0));
    }
    for(size_t i = 1; i < size(); i++) {
      extend_probability_at_site(i, query->get_allele(i));
      if(query->has_span_after(i)) {
        extend_probability_at_span_after(i, query->get_augmentations(i));
      }
    }
  }
  return S;
}

haplotypeMatrix::haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
          haplotypeCohort* cohort, inputHaplotype* query) :
          reference(ref), cohort(cohort), penalties(pen), query(query),
          map(delayMap(cohort->size(), 0)) {
  S = 0;
  R = vector<double>(cohort->size(),0);
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