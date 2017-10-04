#include <cmath>
#include "probability.hpp"

haplotypeMatrix::haplotypeMatrix(const linearReferenceStructure* ref, const penaltySet* pen,
          const haplotypeCohort* cohort) :
          reference(ref), cohort(cohort), penalties(pen),
          map(delayMap(cohort->size(), 0)) {
  S = 0;
  R = vector<double>(cohort->size(), 0);
}

haplotypeMatrix::haplotypeMatrix(const haplotypeMatrix &other, 
            bool copy_map = true) {
	reference = other.reference;
	cohort = other.cohort;
	penalties = other.penalties;
	last_extended = other.last_extended;
	last_span_extended = other.last_span_extended;
	last_allele = other.last_allele;
	S = other.S;
	R = other.R;
	if(copy_map) {
		map = delayMap(other.map);
	} else {
		map = delayMap(cohort->size(), last_extended);
	}
}

haplotypeMatrix::~haplotypeMatrix() {
  
}

void haplotypeMatrix::record_last_extended(alleleValue a) {
  last_extended++;
  last_allele = a;
}

bool haplotypeMatrix::last_extended_is_span() const {
  return (last_extended == last_span_extended);
}

size_t haplotypeMatrix::get_last_site() const {
  return last_extended;
}

delayMap& haplotypeMatrix::get_maps() {
  return map;
}

void haplotypeMatrix::initialize_probability(const inputHaplotype* q) {
  if(q->has_sites()) {
    if(q->has_left_tail()) {
      initialize_probability(q->get_site_index(0), q->get_allele(0),
                q->get_left_tail(), q->get_augmentations(-1));
    } else {
      initialize_probability(q->get_site_index(0), q->get_allele(0));
    }
  } else {
    initialize_probability_at_span(q->get_left_tail(), 
              q->get_augmentations(-1));
  }
}

void haplotypeMatrix::extend_probability_at_site(const inputHaplotype* q, size_t j) {
  extend_probability_at_site(q->get_site_index(j), q->get_allele(j));
}

void haplotypeMatrix::extend_probability_at_span_after(const inputHaplotype* q, 
            size_t j) {
  extend_probability_at_span_after(q->get_site_index(j), 
            q->get_augmentations(j));
}

double haplotypeMatrix::calculate_probability(const inputHaplotype* q) {
  initialize_probability(q);
  if(q->has_span_after(0)) {
    extend_probability_at_span_after(q, 0);
  }
  for(size_t j = 1; j < q->number_of_sites(); j++) {
    extend_probability_at_site(q, j);
    if(q->has_span_after(j)) {
      extend_probability_at_span_after(q, j);
    }
  }
  take_snapshot();
  return S;
}

void haplotypeMatrix::initialize_probability(size_t site_index, alleleValue a,
            size_t left_tail_length, size_t mismatch_count) {
  if(left_tail_length != 0) {
    initialize_probability_at_span(left_tail_length, mismatch_count);
    extend_probability_at_site(site_index, a);
  } else {
    initialize_probability_at_site(site_index, a);
  }
}

void haplotypeMatrix::initialize_probability_at_span(size_t length, 
              size_t mismatch_count) {
  // There is a uniform 1/|H| probability of starting on any given haplotype.
  // All emission probabilities are the same. So all R-values are the same
  length--;
  double common_initial_R = 
            penalties->span_mutation_penalty(length + 1, mismatch_count)
            + penalties->beta(length) - penalties->log_H;
  for(size_t i = 0; i < R.size(); i++) {
    R[i] = common_initial_R;
  }
  S = penalties->span_mutation_penalty(length + 1, mismatch_count)
            + penalties->beta(length);
  
  last_span_extended = -1;
}

void haplotypeMatrix::initialize_probability_at_site(size_t site_index, 
            alleleValue a) {
  vector<size_t> matches = cohort->get_matches(site_index, a);
  vector<size_t> non_matches = 
            cohort->get_non_matches(site_index, a);
  // There are only two possible R-values at this site. There is a uniform
  // 1/|H| probability of starting on any given haplotype; the emission
  // probabilities account for differences in R-value
  double match_initial_value = -penalties->log_H + penalties->one_minus_mu;
  double nonmatch_initial_value = -penalties->log_H + penalties->mu;

  // Set the ones which match the query haplotype
  for(size_t i = 0; i < matches.size(); i++) {
    R[matches[i]] = match_initial_value;
  }
  // And the ones which don't
  for(size_t i = 0; i < non_matches.size(); i++) {
    R[non_matches[i]] = nonmatch_initial_value;
  }

  if(matches.size() == 0) {
    S = penalties->mu;
  } else if(non_matches.size() == 0) {
    S = penalties->one_minus_mu;
  } else {  
    S = -penalties->log_H + 
                logsum(log(matches.size()) + penalties->one_minus_mu,
                       log(non_matches.size()) + penalties->mu);
  }
  record_last_extended(a);
}

void haplotypeMatrix::update_subset_of_Rs(const vector<size_t>& indices,
              bool active_is_match) {
  double correction = penalties->get_minority_map_correction(active_is_match);
  for(size_t i = 0; i < indices.size(); i++) {
    R[indices[i]] = correction + 
                calculate_R(R[indices[i]], map.get_map(indices[i]));
  }
}

void haplotypeMatrix::update_subset_of_Rs(const rowSet& indices,
              bool active_is_match) {
  double correction = penalties->get_minority_map_correction(active_is_match);
  for(size_t i = 0; i < indices.size(); i++) {
    R[indices[i]] = correction + 
                calculate_R(R[indices[i]], map.get_map(indices[i]));
  }
}

void haplotypeMatrix::fast_update_S(const vector<size_t>& indices,
              bool active_is_match) {
  vector<double> summands = vector<double>(indices.size(), 0);
  for(size_t i = 0; i < indices.size(); i++) {
    summands[i] = R[indices[i]];
  }
  penalties->update_S(S, summands, active_is_match);
}

void haplotypeMatrix::fast_update_S(const rowSet& indices,
              bool active_is_match) {
  vector<double> summands = vector<double>(indices.size(), 0);
  for(size_t i = 0; i < indices.size(); i++) {
    summands[i] = R[indices[i]];
  }
  penalties->update_S(S, summands, active_is_match);
}

void haplotypeMatrix::extend_probability_at_site(size_t site_index,
            alleleValue a) {
  bool match_is_rare = cohort->match_is_rare(site_index, a);
  DPUpdateMap current_map = penalties->get_current_map(S, match_is_rare);
  vector<size_t> active_rows = cohort->get_active_rows(site_index, a);
  extend_probability_at_site(current_map, active_rows, match_is_rare, a);
}

void haplotypeMatrix::extend_probability_at_span_after(size_t site_index,
            size_t mismatch_count = 0) {
  size_t length = reference->span_length_after(site_index);
  extend_probability_at_span_after_anonymous(length, mismatch_count);
}

void haplotypeMatrix::take_snapshot() {
  if(last_extended >= 0) {
    // step forward all delayMap slots which are not up to date
    map.hard_update_all();
    size_t j = last_extended;
    bool reference_is_homogenous = 
              (cohort->number_matching(j, last_allele) == 0 ||
              cohort->number_not_matching(j, last_allele) == 0);
    if(reference_is_homogenous || last_extended_is_span()) {
      for(size_t i = 0; i < cohort->size(); i++) {
        R[i] = calculate_R(R[i], map.get_map(i));
      }
    } else if(cohort->match_is_rare(j, last_allele)) {
      // we have already calculated R at all matching rows
      vector<size_t> non_matches = cohort->get_non_matches(j, last_allele);
      for(size_t i = 0; i < non_matches.size(); i++) {
        R[non_matches[i]] =
                  calculate_R(R[non_matches[i]], map.get_map(non_matches[i]));
      }
    } else {
      // we have already calculated R at all non-matching rows
      vector<size_t> matches = cohort->get_matches(j, last_allele);
      for(size_t i = 0; i < matches.size(); i++) {
        R[matches[i]] =
                  calculate_R(R[matches[i]], map.get_map(matches[i]));
      }
    }
    // since all R-values are up to date, we do not need entries in the delayMap
    // therefore we can clear them all and replace them with the identity map
    map.hard_clear_all();
  }
} 

double haplotypeMatrix::prefix_likelihood() const {
  return S;
}

double haplotypeMatrix::partial_likelihood_by_row(size_t row) const {
  return R[row];
}

double calculate_R(double oldR, const DPUpdateMap& map) {
  return map.of(oldR);
}

double calculate_R(double oldR, double coefficient, double constant) {
  return calculate_R(oldR, DPUpdateMap(coefficient, constant));
}

void haplotypeMatrix::extend_probability_at_site(const DPUpdateMap& current_map, 
            const vector<size_t>& active_rows, bool match_is_rare, 
            alleleValue a) {
  map.add_map_for_site(current_map);
  if(active_rows.size() == 0 && match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->mu + S + penalties->beta_value;
  } else if(active_rows.size() == 0 && !match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->one_minus_mu + S + penalties->beta_value;
  } else {
    map.update_map_with_active_rows(active_rows);
    update_subset_of_Rs(active_rows, match_is_rare);
    fast_update_S(active_rows, match_is_rare);
    map.reset_rows(active_rows);
  }
  record_last_extended(a);
  return;
}

void haplotypeMatrix::extend_probability_at_site(const DPUpdateMap& current_map, 
            const rowSet& active_rows, bool match_is_rare, 
            alleleValue a) {
  map.add_map_for_site(current_map);
  if(active_rows.size() == 0 && match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->mu + S + penalties->beta_value;
  } else if(active_rows.size() == 0 && !match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->one_minus_mu + S + penalties->beta_value;
  } else {
    map.update_map_with_active_rows(active_rows);
    update_subset_of_Rs(active_rows, match_is_rare);
    fast_update_S(active_rows, match_is_rare);
    map.reset_rows(active_rows);
  }
  record_last_extended(a);
  return;
}

void haplotypeMatrix::extend_probability_at_site(
            const vector<size_t>& active_rows, bool match_is_rare, 
            alleleValue a) {
  DPUpdateMap current_map = penalties->get_current_map(S, match_is_rare);
  extend_probability_at_site(current_map, active_rows, match_is_rare, a);
}

void haplotypeMatrix::extend_probability_at_site(
            const rowSet& active_rows, bool match_is_rare, 
            alleleValue a) {
  DPUpdateMap current_map = penalties->get_current_map(S, match_is_rare);
  extend_probability_at_site(current_map, active_rows, match_is_rare, a);
}

void haplotypeMatrix::extend_probability_at_span_after_anonymous(size_t l, 
            size_t mismatch_count) {
  double m = penalties->span_mutation_penalty(l, mismatch_count);
  map.update_map_with_span(m + penalties->alpha(l), 
              S + penalties->span_coefficient(l) - penalties->alpha(l));
  S = m + S + penalties->beta(l);
  last_span_extended = last_extended;
}