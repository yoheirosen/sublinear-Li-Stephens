#include <cmath>
#include "math.hpp"
#include "probability.hpp"
#include <iostream>

#ifdef TIME_PROBABILITY_INTERNALS
#include <sys/time.h>
#endif

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

#ifdef TIME_PROBABILITY_INTERNALS
void fastFwdAlgState::set_timers(double* total, double* readwrite, double* delay) {
  t_total = total;
  t_readwrite = readwrite;
  t_delay = delay;
  *t_total = 0;
  *t_readwrite = 0;
  *t_delay = 0;
}
#endif

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
      extend_probability_at_site(q->get_site_index(0), q->get_allele(0));
    } else {
      initialize_probability_at_site(q->get_site_index(0), q->get_allele(0));
    }
  } else {
    initialize_probability_at_span(q->get_left_tail());
  }
}

double fastFwdAlgState::calculate_probability(const inputHaplotype* observed) {
#ifdef DEBUG
  vector<vector<double> > fast_values(observed->number_of_sites(), vector<double>(cohort->get_n_haplotypes(), 1);
#endif

#ifdef TIME_PROBABILITY_INTERNALS
  struct timeval timer1, timer2;
  gettimeofday(&timer1, NULL);
#endif
  
  initialize_probability(observed);
  
#ifdef TIME_PROBABILITY_INTERNALS
  gettimeofday(&timer2, NULL);
  *t_total += (double) (timer2.tv_usec - timer1.tv_usec) / 1000000 + (double) (timer2.tv_sec - timer1.tv_sec);
#endif

#ifdef DEBUG
  rowSet& indices = reference->get_active_rowSet()
  rowSet::const_iterator rows_begin = indices.begin();
  rowSet::const_iterator rows_end = indices.end();
  for(rowSet::const_iterator it = rows_begin; it != rows_end; ++it) {
    size_t row = *it;
    fast_values[0][row] = rows[row];
  }
#endif

  if(observed->has_span_after(0)) {
    extend_probability_at_span_after_abstract(observed->get_span_after(0));
  }
  
  for(size_t i = 1; i < observed->number_of_sites(); i++) {
    extend_probability_at_site(observed->get_site_index(i), observed->get_allele(i));
    
#ifdef TIME_PROBABILITY_INTERNALS
    gettimeofday(&timer2, NULL);
    *t_total += (double) (timer2.tv_usec - timer1.tv_usec) / 1000000 + (double) (timer2.tv_sec - timer1.tv_sec);
#endif
#ifdef DEBUG
    for(rowSet::const_iterator it = rows_begin; it != rows_end; ++it) {
      size_t row = *it;
      fast_values[i][row] = rows[row];
    }
#endif
#ifdef TIME_PROBABILITY_INTERNALS
      gettimeofday(&timer1, NULL);
#endif
    
    if(observed->has_span_after(i)) {
      extend_probability_at_span_after_abstract(observed->get_span_after(i));
    }
  }

#ifdef DEBUG
  for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
    for(size_t i = 0; i < observed->number_of_sites(); i++) {
      cerr << "[";
      if(fast_values[i][j] > 0) { cerr << "*"; } else { cerr << fast_values[i][j]; }
      cerr << "]\t";
    }
    cerr << endl;
  }
#endif
  
#ifdef TIME_PROBABILITY_INTERNALS
  gettimeofday(&timer2, NULL);
  *t_total += (double) (timer2.tv_usec - timer1.tv_usec) / 1000000 + (double) (timer2.tv_sec - timer1.tv_sec);
  *t_readwrite = *t_readwrite - *t_delay;
#endif

  return sum;
}

void fastFwdAlgState::initialize_probability_at_span(size_t length) {
  // There is a uniform 1/|k| probability of starting on any given haplotype.
  // All emission probabilities are the same. So all row-values are the same
  double common_initial_row = penalties->span_mutation_penalty(length) - penalties->log_k;
  std::fill(rows.begin(), rows.end(), common_initial_row);
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
                logmath::logsum(log(cohort->number_matching(site_index, a)) + penalties->one_minus_mu,
                       log(cohort->number_not_matching(site_index, a)) + penalties->mu);
  }
  record_last_extended(a);
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
  
#ifdef TIME_PROBABILITY_INTERNALS  
  struct timeval timer1, timer2;
  gettimeofday(&timer1, NULL);
  map.update_evaluate_and_move_rows(active_row_indices, rows, penalties->minority_correction(match_is_rare), t_delay);
#else 

  map.update_evaluate_and_move_rows(active_row_indices, rows, penalties->minority_correction(match_is_rare));

#endif
#ifdef TIME_PROBABILITY_INTERNALS  
  gettimeofday(&timer2, NULL);
  *t_readwrite += (double) (timer2.tv_usec - timer1.tv_usec) / 1000000 + (double) (timer2.tv_sec - timer1.tv_sec);
#endif

  update_sum(active_row_indices, match_is_rare);
  record_last_extended(a);
  return;
}

slowFwdSolver::slowFwdSolver(siteIndex* ref, const penaltySet* pen, const haplotypeCohort* haplotypes) :
            reference(ref), penalties(pen), cohort(haplotypes) {
}

double slowFwdSolver::calculate_probability_quadratic(const inputHaplotype* q) {
  const vector<alleleValue>& alleles = q->get_alleles();
  size_t start_site = q->get_start_site();
  sum = q->has_left_tail() ? penalties->span_mutation_penalty(q->get_left_tail()) : 0;
  double initial_value = -penalties->log_k + sum;
  rows = vector<double>(cohort->get_n_haplotypes(), initial_value);
  for(size_t j = 0; j < rows.size(); j++) {
    bool matches = (cohort->allele_at(start_site, j) == alleles[0]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    rows[j] += emission;
  }
  sum = logmath::log_big_sum(rows);

  vector<double> last_rows;
  double same_transition = log1p(-exp(penalties->rho) * (penalties->k - 1));
  for(size_t i = 1; i < q->number_of_sites(); i++) {
    // do span stuff
    if(q->has_span_after(i - 1)) {
      size_t l = q->get_span_after(i - 1);
      double coefficient = penalties->pow_rho_c(l);
      double constant = penalties->span_polynomial(l) + sum;
      for(size_t j = 0; j < rows.size(); j++) {
        rows[j] = penalties->span_mutation_penalty(l) + logmath::logsum(coefficient + rows[j], constant);
      }
    }
    sum = logmath::log_big_sum(rows);
    last_rows = rows;
    for(size_t j = 0; j < rows.size(); j++) {
      vector<double> temp = last_rows;
      for(size_t k = 0; k < temp.size(); k++) {
        temp[k] += j == k ? same_transition : penalties->rho;
      }
      rows[j] = logmath::log_big_sum(temp);
      bool matches = (cohort->allele_at(start_site + i, j) == alleles[i]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      rows[j] += emission;
    }
    sum = logmath::log_big_sum(rows);
  }
  if(q->has_span_after(q->number_of_sites() - 1)) {
    size_t l = q->get_span_after(q->number_of_sites() - 1);
    double coefficient = penalties->pow_rho_c(l);
    double constant = penalties->span_polynomial(l) + sum;
    for(size_t j = 0; j < rows.size(); j++) {
      rows[j] = penalties->span_mutation_penalty(l) + logmath::logsum(coefficient + rows[j], constant);
    }
  }
  
  sum = logmath::log_big_sum(rows);
  return sum;
}

double slowFwdSolver::calculate_probability_linear(const inputHaplotype* q) {
  initialize_linear(q);
  
#ifdef DEBUG
  vector<vector<double> > slow_values(q->number_of_sites(), vector<double>(cohort->get_n_haplotypes()));
  for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
    slow_values[0][j] = rows[j];
  }
#endif
  
  for(size_t i = 1; i < q->number_of_sites(); i++) {
    if(q->has_span_after(i - 1)) {
      extend_span_linear(q, i - 1);
    }
    extend_site_linear(q, i);
    
#ifdef DEBUG
    for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
      slow_values[i][j] = rows[j];
    }
#endif

  }
  if(q->has_span_after(q->number_of_sites() - 1)) {
    extend_span_linear(q, q->number_of_sites() - 1);
  }
  
#ifdef DEBUG
  for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
    for(size_t i = 0; i < q->number_of_sites(); i++) {
      cerr << "[" << slow_values[i][j] << "]\t";
    }
    cerr << endl;
  }
#endif

  return sum;
}

void slowFwdSolver::initialize_linear(const inputHaplotype* q) {
  size_t start_site = q->get_start_site();
  sum = q->has_left_tail() ? penalties->span_mutation_penalty(q->get_left_tail()) : 0;
  rows = vector<double>(cohort->get_n_haplotypes(), -penalties->log_k + sum);
  for(size_t j = 0; j < rows.size(); j++) {
    bool matches = (cohort->allele_at(start_site, j) == q->get_allele(0));
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    rows[j] += emission;
  }
  sum = logmath::log_big_sum(rows);
}

void slowFwdSolver::extend_site_linear(const inputHaplotype* q, size_t site) {
  size_t start_site = q->get_start_site();
  for(size_t j = 0; j < rows.size(); j++) {
    rows[j] = logmath::logsum(sum + penalties->rho, rows[j] + penalties->row_coefficient); 
    bool matches = (cohort->allele_at(start_site + site, j) == q->get_allele(site));
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    rows[j] += emission;
  }
  // sum = logmath::log_big_sum(rows);
  if(!cohort->match_is_rare(start_site + site, q->get_allele(site))) {
    if(cohort->number_not_matching(start_site + site, q->get_allele(site)) != 0) {
      const rowSet& indices = cohort->get_active_rowSet(start_site + site, q->get_allele(site));
      rowSet::const_iterator it = indices.begin();
      rowSet::const_iterator rows_end = indices.end();
      // vector<double> to_correct;
      // for(it; it != rows_end; ++it) {
      //   size_t row = *it;
      //   to_correct.push_back(rows[row]);
      // }
      sum = logmath::logdiff(penalties->one_minus_mu + sum, penalties->one_minus_2mu - penalties->mu + logmath::log_big_sum(it, rows_end, rows));
    } else {
      sum += penalties->one_minus_mu;
    }
  } else {
    if(cohort->number_matching(start_site + site, q->get_allele(site)) != 0) {
      const rowSet& indices = cohort->get_active_rowSet(start_site + site, q->get_allele(site));
      rowSet::const_iterator it = indices.begin();
      rowSet::const_iterator rows_end = indices.end();
      // vector<double> to_correct;
      // for(it; it != rows_end; ++it) {
      //   size_t row = *it;
      //   to_correct.push_back(rows[row]);
      // }
      sum = logmath::logsum(penalties->mu + sum, penalties->one_minus_2mu - penalties->one_minus_mu + logmath::log_big_sum(it, rows_end, rows));
    } else {
      sum += penalties->mu;
    }
  }
}

void slowFwdSolver::extend_span_linear(const inputHaplotype* q, size_t site) {
  if(q->has_span_after(site)) {
    size_t l = q->get_span_after(site);
    double m = penalties->span_mutation_penalty(l);
    double coefficient = penalties->pow_rho_c(l);
    double constant = penalties->span_polynomial(l) + sum;
    for(size_t j = 0; j < rows.size(); j++) {
      rows[j] = penalties->span_mutation_penalty(l) + logmath::logsum(coefficient + rows[j], constant);
    }
    sum += penalties->span_mutation_penalty(l);
    // sum = logmath::log_big_sum(rows);
  }
}