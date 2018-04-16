#include <cmath>
#include "probability.hpp"
#include <iostream>

const fastFwdAlgState::dp_column_t fastFwdAlgState::SITE_UNEXTENDED = -1;
const fastFwdAlgState::dp_column_t fastFwdAlgState::SPAN_UNEXTENDED = -2;
const fastFwdAlgState::dp_column_t fastFwdAlgState::INITIAL_SPAN = -1;

fastFwdAlgState::fastFwdAlgState(siteIndex* reference, const penaltySet* penalties, const haplotypeCohort* cohort) :
          reference(reference), cohort(cohort), penalties(penalties), map(delayedEvalMap(cohort->get_n_haplotypes())) {
  sum = 0;
  rows = vector<double>(cohort->get_n_haplotypes(), 0);
}

fastFwdAlgState::fastFwdAlgState(const fastFwdAlgState &other, bool copy_map = true) {
	reference = other.reference;
	cohort = other.cohort;
	penalties = other.penalties;
	last_extended = other.last_extended;
	last_span_extended = other.last_span_extended;
	last_allele = other.last_allele;
	sum = other.sum;
	rows = other.rows;
	if(copy_map) {
		map = delayedEvalMap(other.map);
	} else {
		map = delayedEvalMap(cohort->get_n_haplotypes());
	}
}

fastFwdAlgState::~fastFwdAlgState() {
  
}

inline void fastFwdAlgState::record_last_extended(alleleValue a) {
  last_extended++;
  last_allele = a;
}

inline bool fastFwdAlgState::last_extended_is_span() const {
  return (last_extended == last_span_extended);
}

inline size_t fastFwdAlgState::get_last_site() const {
  return last_extended;
}

void fastFwdAlgState::initialize_probability(const inputHaplotype* q) {
  if(q->has_sites()) {
    if(q->has_left_tail()) {
      initialize_probability_at_span(q->get_left_tail());
    } else {
      initialize_probability_at_site(q->get_site_index(0), q->get_allele(0));
    }
  } else {
    initialize_probability_at_span(q->get_left_tail());
  }
}

double fastFwdAlgState::calculate_probability(const inputHaplotype* observed) {
  initialize_probability(observed);
  if(observed->has_span_after(0)) {
    extend_probability_at_span_after_abstract(observed->get_span_after(0));
  }
  for(size_t i = 1; i < observed->number_of_sites(); i++) {
    extend_probability_at_site(observed->get_site_index(i), observed->get_allele(i));
    if(observed->has_span_after(i)) {
      extend_probability_at_span_after_abstract(observed->get_span_after(i));
    }
  }
  return sum;
}

void fastFwdAlgState::initialize_probability_at_span(size_t length) {
  // There is a uniform 1/|k| probability of starting on any given haplotype.
  // All emission probabilities are the same. So all row-values are the same
  double common_initial_row = penalties->span_mutation_penalty(length) - penalties->log_k;
  for(size_t i = 0; i < rows.size(); i++) {
    rows[i] = common_initial_row;
  }
  sum = penalties->span_mutation_penalty(length);  
  last_span_extended = -1;
}

void fastFwdAlgState::initialize_probability_at_site(size_t site_index, alleleValue a) {
  // There are only two possible row-values at this site. There is a uniform
  // 1/|k| probability of starting on any given haplotype; the emission
  // probabilities account for differences in row-value
  double match_initial_value = -penalties->log_k + penalties->one_minus_mu;
  double nonmatch_initial_value = -penalties->log_k + penalties->mu;

  double active_value = cohort->match_is_rare(site_index, a) ? match_initial_value : nonmatch_initial_value;
  double default_value = cohort->match_is_rare(site_index, a) ? nonmatch_initial_value : match_initial_value;
  
  std::fill(rows.begin(), rows.end(), default_value);
  
  if(cohort->number_active(site_index, a) != 0) {
    const rowSet& active_row_indices = cohort->get_active_rowSet(site_index, a);
    rowSet::const_iterator it = active_row_indices.begin();
    rowSet::const_iterator rows_end = active_row_indices.end();
    for(it; it != rows_end; ++it) {
      rows[*it] = active_value;
    }
  }

  if(cohort->number_matching(site_index, a) == 0) {
    sum = penalties->mu;
  } else if(cohort->number_not_matching(site_index, a) == 0) {
    sum = penalties->one_minus_mu;
  } else {  
    sum = -penalties->log_k + 
                logsum(log(cohort->number_matching(site_index, a)) + penalties->one_minus_mu,
                       log(cohort->number_not_matching(site_index, a)) + penalties->mu);
  }
  record_last_extended(a);
}

void fastFwdAlgState::update_subset_of_rows(const rowSet& indices,
              bool active_is_match) {
  double correction = penalties->minority_correction(active_is_match);
  map.update_evaluate_and_move_rows(indices, rows, correction);
}

void fastFwdAlgState::update_sum(const rowSet& indices, bool active_is_match) {
  penalties->update_sum(sum, rows, indices.begin(), indices.end(), active_is_match);
}

void fastFwdAlgState::extend_probability_at_site(size_t site_index, alleleValue a) {
  bool match_is_rare = cohort->match_is_rare(site_index, a);
  DPUpdateMap current_majority_map = penalties->current_map(sum, match_is_rare);
  const rowSet& active_row_indices = cohort->get_active_rowSet(site_index, a);
  extend_probability_at_site(current_majority_map, active_row_indices, match_is_rare, a);
}

void fastFwdAlgState::take_snapshot() {
  if(last_extended >= 0) {
    map.hard_update_all();
    size_t j = last_extended;
    bool reference_is_homogeneous = (cohort->number_matching(j, last_allele) == 0 || cohort->number_not_matching(j, last_allele) == 0);
    if(reference_is_homogeneous || last_extended_is_span()) {
      for(size_t i = 0; i < cohort->get_n_haplotypes(); i++) {
        rows[i] = calculate_row(rows[i], map.get_map(i));
      }
    } else {
      vector<bool> already_calculated(cohort->get_n_haplotypes(), false);
      rowSet::const_iterator it = cohort->get_active_rowSet(j, last_allele).begin();
      rowSet::const_iterator end = cohort->get_active_rowSet(j, last_allele).end();
      for(it; it != end; ++it) {
        already_calculated[*it] = true;
      }
      for(size_t i = 0; i < cohort->get_n_haplotypes(); i++) {
        if(!already_calculated[i]) {
          rows[i] = calculate_row(rows[i], map.get_map(i));
        }
      }
    }
  }
} 

double calculate_row(double old_row, const DPUpdateMap& map) {
  return map.of(old_row);
}

double calculate_row(double old_row, double coefficient, double constant) {
  return calculate_row(old_row, DPUpdateMap(coefficient, constant));
}

void fastFwdAlgState::extend_probability_at_span_after_abstract(size_t l) {
  double m = penalties->span_mutation_penalty(l);
  double coefficient = penalties->pow_rho_c(l);
  double constant = penalties->span_polynomial(l) + sum;
  map.extend_value_only(DPUpdateMap(m + coefficient, constant - coefficient));
  sum += m;
  last_span_extended = last_extended;
}

void fastFwdAlgState::extend_probability_at_site(const DPUpdateMap& current_majority_map, 
                                                 const rowSet& active_row_indices, 
                                                 bool match_is_rare, 
                                                 alleleValue a) {
  if(active_row_indices.empty()) {
    if(match_is_rare) {
      // separate case to avoid log-summing "log 0"
      sum += penalties->mu;
    } else {
      // separate case to avoid log-summing "log 0"
      sum += penalties->one_minus_mu;
    }
    map.extend_value_only(current_majority_map);
    return;
  }

  map.extend_by_new_step(current_majority_map);
  // update rows by composing with everything *including this map*
  // don't update least recently updated
  // then use same loop to evaluate as well as move
  map.update_evaluate_and_move_rows(active_row_indices, rows, penalties->minority_correction(match_is_rare));
  update_sum(active_row_indices, match_is_rare);
  record_last_extended(a);
  return;
}

slowFwdSolver::slowFwdSolver(siteIndex* ref, const penaltySet* pen, const haplotypeCohort* haplotypes) :
            reference(ref), penalties(pen), cohort(haplotypes) {
}
            
double slowFwdSolver::calculate_probability_quadratic(const vector<alleleValue>& q, size_t start_site = 0) {
  rows = vector<double>(cohort->get_n_haplotypes(), -penalties->log_k);
  for(size_t j = 0; j < rows.size(); j++) {
    bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    rows[j] += emission;
  }
  vector<double> last_rows;
  double same_transition = log1p(-exp(penalties->rho) * (penalties->k - 1));
  for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
    last_rows = rows;
    if(reference->has_span_after(i - 1)) {
      size_t l = reference->span_length_after(i - 1);
      double coefficient = penalties->span_mutation_penalty(l) + penalties->composed_row_coefficient(l);
      double constant = penalties->span_coefficient(l) - penalties->composed_row_coefficient(l);
      for(size_t j = 0; j < rows.size(); j++) {
        rows[j] = coefficient + logsum(rows[j], constant);
      }
    }
    for(size_t j = 0; j < rows.size(); j++) {
      vector<double> temp = last_rows;
      for(size_t k = 0; k < temp.size(); k++) {
        temp[k] += j == k ? same_transition : penalties->rho;
      }
      rows[j] = log_big_sum(temp);
      bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      rows[j] += emission;
    }
  }
  sum = log_big_sum(rows);
  return sum;
}

double slowFwdSolver::calculate_probability_linear(const vector<alleleValue>& q, size_t start_site = 0) {
  rows = vector<double>(cohort->get_n_haplotypes(), -penalties->log_k);
  for(size_t j = 0; j < rows.size(); j++) {
    bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    rows[j] += emission;
  }
  vector<double> last_rows;
  for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
    last_rows = rows;
    sum = log_big_sum(last_rows);
    // do span stuff
    if(reference->has_span_after(i - 1)) {
      size_t l = reference->span_length_after(i - 1);
      double coefficient = penalties->span_mutation_penalty(l) + penalties->composed_row_coefficient(l);
      double constant = penalties->span_coefficient(l) + sum - penalties->composed_row_coefficient(l);
      for(size_t j = 0; j < rows.size(); j++) {
        rows[j] = coefficient + logsum(rows[j], constant);
      }
    }
    for(size_t j = 0; j < rows.size(); j++) {
      rows[j] = logsum(sum + penalties->rho, last_rows[j] + penalties->row_coefficient); 
      bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      rows[j] += emission;
    }
  }
  sum = log_big_sum(rows);
  return sum;
}

double slowFwdSolver::calculate_probability_quadratic(const inputHaplotype* observed_haplotype) {
  return calculate_probability_quadratic(observed_haplotype->get_alleles(), observed_haplotype->get_start_site());
}

double slowFwdSolver::calculate_probability_linear(const inputHaplotype* observed_haplotype) {
  return calculate_probability_linear(observed_haplotype->get_alleles(), observed_haplotype->get_start_site());
}

// slowFwdSolver::sequence_statistic_array slowFwdSolver::sequence_statistics(const vector<alleleValue>& q, size_t start_site) {
//   for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
//     alleleValue test = q.at(i - start_site);
//     if(test > 5) {
//       cout << i << " " << i - start_site << " " << test << " " << sizeof(alleleValue) << endl;
//       throw runtime_error("invalid allele in observed haplotype");
//     }
//   }
//   
//   vector<size_t> unique_values(q.size(), 1);
//   vector<size_t> unique_active_values(q.size(), 0);
//   vector<size_t> active_row_indices(q.size(), 0);
//   vector<size_t> history_length(q.size(), 0);
//   vector<size_t> active_history_length(q.size(), 0);
//   vector<size_t> map_classes(q.size(), 0);
//   vector<size_t> active_map_classes(q.size(), 0);
//   
//   vector<size_t> times_active(cohort->get_n_haplotypes(), 
//   0);
//   vector<size_t> last_updated(cohort->get_n_haplotypes(), start_site);
//   rows = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H);
//   for(size_t j = 0; j < rows.size(); j++) {
//     bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
//     double emission = matches ? penalties->one_minus_mu : penalties->mu;
//     rows[j] += emission;
//   }
//   vector<double> last_rows;
//   for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
//     last_rows = rows;
//     sum = log_big_sum(last_rows);
//     // do span stuff
//     if(reference->has_span_after(i - 1)) {
//       size_t l = reference->span_length_after(i - 1);
//       double coefficient = penalties->span_mutation_penalty(l, 0) + penalties->composed_row_coefficient(l);
//       double constant = penalties->span_coefficient(l) - penalties->composed_row_coefficient(l);
//       for(size_t j = 0; j < rows.size(); j++) {
//         rows[j] = coefficient + logsum(rows[j], constant);
//       }
//     }
//     for(size_t j = 0; j < rows.size(); j++) {
//       rows[j] = logsum(sum + penalties->rho, last_rows[j] + penalties->row_coefficient); 
//       bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
//       double emission = matches ? penalties->one_minus_mu : penalties->mu;
//       rows[j] += emission;
//     }
//     sum = log_big_sum(rows);
//     // DATA gathering  
// 
//     double eps = 0.000000001;
//     // VALUES --------------------------------------------------------------
//     for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(fabs(rows.at(j) - rows.at(k)) < eps) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++unique_values.at(i - start_site);
//       }
//     }
//     vector<size_t> i_active_row_indices = cohort->get_active_row_indices(i, q.at(i - start_site));
//     for(size_t j = 0; j < i_active_row_indices.size(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(fabs(rows.at(i_active_row_indices.at(j)) - rows.at(i_active_row_indices.at(k))) < eps) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++unique_active_values.at(i - start_site);
//       }
//     }
//     
//     // ACTIVE ROWS --------------------------------------------------------------
//     active_row_indices.at(i - start_site) = i_active_row_indices.size();
//     // HISTORY LENGTH -----------------------------------------------------------
//     size_t oldest_active = i;
//     size_t n_active_mapclasses = 0;
//     size_t oldest = i;
//     size_t n_mapclasses = 0;
//     for(size_t j = 0; j < i_active_row_indices.size(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(last_updated.at(i_active_row_indices.at(j)) == last_updated.at(i_active_row_indices.at(k))) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++n_active_mapclasses;
//       }
//       if(last_updated.at(i_active_row_indices.at(j)) < oldest_active) {
//         oldest_active = last_updated.at(i_active_row_indices.at(j));
//       }
//     }
//     
//     for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(last_updated.at(j) == last_updated.at(k)) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++n_mapclasses;
//       }
//       if(last_updated.at(j) < oldest) {
//         oldest = last_updated.at(j);
//       }
//     }
//     
//     history_length.at(i - start_site) = i - oldest;
//     active_history_length.at(i -start_site) = i - oldest_active;
//     map_classes.at(i - start_site) = n_mapclasses;
//     active_map_classes.at(i - start_site) = n_active_mapclasses;
//     
//     // Track activity
//     for(size_t j = 0; j < i_active_row_indices.size(); j++) {
//       last_updated.at(i_active_row_indices.at(j)) = i;
//       ++times_active.at(i_active_row_indices.at(j));
//     }
//   }
//   
//   size_t max_u_value = 1;
//   double avg_u_value = 0;
//   size_t max_u_active_value = 1;
//   double avg_u_active_value = 0;
//   size_t max_active_row_indices = 0;
//   double avg_active_row_indices = 0;
//   size_t max_history_length = 0;
//   double avg_history_length = 0;
//   size_t max_active_history_length = 0;
//   double avg_active_history_length = 0;
//   size_t max_map_classes = 0;
//   double avg_map_classes = 0;
//   size_t max_active_map_classes = 0;
//   double avg_active_map_classes = 0;
//   size_t max_times_active = 0;
//   double avg_times_active = 0;
//     
//   for(size_t i = 0; i < unique_values.size(); i++) {
//     if(unique_values.at(i) > max_u_value) {
//       max_u_value = unique_values.at(i);
//     }
//     if(unique_active_values.at(i) > max_u_active_value) {
//       max_u_active_value = unique_active_values.at(i);
//     }
//     if(active_row_indices.at(i) > max_active_row_indices) {
//       max_active_row_indices = active_row_indices.at(i);
//     }
//     if(history_length.at(i) > max_history_length) {
//       max_history_length = history_length.at(i);
//     }
//     if(active_history_length.at(i) > max_active_history_length) {
//       max_active_history_length = active_history_length.at(i);
//     }
//     if(active_map_classes.at(i) > max_active_map_classes) {
//       max_active_map_classes = active_map_classes.at(i);
//     }
//     if(map_classes.at(i) > max_map_classes) {
//       max_map_classes = map_classes.at(i);
//     }
//     
//     avg_u_value += unique_values.at(i);
//     avg_u_active_value += unique_active_values.at(i);
//     avg_active_row_indices += active_row_indices.at(i);
//     avg_history_length += history_length.at(i);
//     avg_active_history_length += active_history_length.at(i);
//     avg_active_map_classes += active_map_classes.at(i);
//     avg_map_classes += map_classes.at(i);
//   }
//   for(size_t j = 0; j < times_active.size(); j++) {
//     if(times_active.at(j) > max_times_active) {
//       max_times_active = times_active.at(j);
//     }
//     avg_times_active += times_active.at(j);
//   }
//   
//   size_t common_denom = unique_values.size() - 1;
//   avg_u_value = avg_u_value/common_denom;
//   avg_u_active_value = avg_u_active_value/common_denom;
//   avg_active_row_indices = avg_active_row_indices/common_denom;
//   avg_history_length = avg_history_length/common_denom;
//   avg_active_history_length = avg_active_history_length/common_denom;
//   avg_map_classes = avg_map_classes/common_denom;
//   avg_active_map_classes = avg_active_map_classes/common_denom;
//   avg_times_active = avg_times_active/times_active.size();
//   
//   vector<size_t> to_return_max = {
//     max_u_value,
//     max_u_active_value,
//     max_active_row_indices,
//     max_history_length,
//     max_active_history_length,
//     max_map_classes,
//     max_active_map_classes,
//     max_times_active
//   };
// 
//   vector<double> to_return_avg = {
//     avg_u_value,
//     avg_u_active_value,
//     avg_active_row_indices,
//     avg_history_length,
//     avg_active_history_length,
//     avg_map_classes,
//     avg_active_map_classes,
//     avg_times_active
//   };
//   
//   return make_pair(to_return_avg, to_return_max);
// }