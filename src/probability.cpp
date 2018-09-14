#include <cmath>
#include "probability.hpp"

struct liStephensModel{
  liStephensModel(siteIndex* reference, haplotypeCohort* cohort, const penaltySet* penalties);
  siteIndex* reference;
  haplotypeCohort* cohort;
  const penaltySet* penalties;
};

fastFwdAlgState::fastFwdAlgState(siteIndex* reference, const penaltySet* penalties, const haplotypeCohort* cohort) :
          reference(reference), cohort(cohort), penalties(penalties), map(lazyEvalMap(cohort->get_n_haplotypes(), 0)) {
  S = 0;
  R = vector<double>(cohort->get_n_haplotypes(), 0);
}

fastFwdAlgState::fastFwdAlgState(const fastFwdAlgState &other, bool copy_map) {
	reference = other.reference;
	cohort = other.cohort;
	penalties = other.penalties;
	last_extended = other.last_extended;
	last_span_extended = other.last_span_extended;
	last_allele = other.last_allele;
	S = other.S;
	R = other.R;
	if(copy_map) {
		map = lazyEvalMap(other.map);
	} else {
		map = lazyEvalMap(cohort->get_n_haplotypes(), last_extended);
	}
}

fastFwdAlgState::~fastFwdAlgState() {
  
}

void fastFwdAlgState::record_last_extended(alleleValue a) {
  last_extended++;
  last_allele = a;
}

bool fastFwdAlgState::last_extended_is_span() const {
  return (last_extended == last_span_extended);
}

size_t fastFwdAlgState::get_last_site() const {
  return last_extended;
}

lazyEvalMap& fastFwdAlgState::get_maps() {
  return map;
}

void fastFwdAlgState::initialize_probability(const inputHaplotype* q) {
  if(q->has_sites()) {
    if(q->has_left_tail()) {
      initialize_probability(q->get_site_index(0), q->get_allele(0),
                q->get_left_tail(), q->get_n_novel_SNVs(-1));
    } else {
      initialize_probability(q->get_site_index(0), q->get_allele(0));
    }
  } else {
    initialize_probability_at_span(q->get_left_tail(), 
              q->get_n_novel_SNVs(-1));
  }
}

void fastFwdAlgState::extend_probability_at_site(const inputHaplotype* q, size_t j) {
  extend_probability_at_site(q->get_site_index(j), q->get_allele(j));
}

void fastFwdAlgState::extend_probability_at_span_after(const inputHaplotype* q, 
            size_t j) {
  extend_probability_at_span_after_anonymous(q->get_span_after(j), q->get_n_novel_SNVs(j));
}

double fastFwdAlgState::calculate_probability(const inputHaplotype* q) {
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
  return S;
}

void fastFwdAlgState::initialize_probability(size_t site_index, alleleValue a,
            size_t left_tail_length, size_t mismatch_count) {
  if(left_tail_length != 0) {
    initialize_probability_at_span(left_tail_length, mismatch_count);
    extend_probability_at_site(site_index, a);
  } else {
    initialize_probability_at_site(site_index, a);
  }
}

void fastFwdAlgState::initialize_probability_at_span(size_t length, 
              size_t mismatch_count) {
  // There is a uniform 1/|H| probability of starting on any given haplotype.
  // All emission probabilities are the same. So all R-values are the same
  double common_initial_R = penalties->span_mutation_penalty(length, mismatch_count) - penalties->log_H;
  for(size_t i = 0; i < R.size(); i++) {
    R[i] = common_initial_R;
  }
  S = penalties->span_mutation_penalty(length, mismatch_count);  
  last_span_extended = -1;
}

void fastFwdAlgState::initialize_probability_at_site(size_t site_index, 
            alleleValue a) {
  // There are only two possible R-values at this site. There is a uniform
  // 1/|H| probability of starting on any given haplotype; the emission
  // probabilities account for differences in R-value
  double match_initial_value = -penalties->log_H + penalties->one_minus_mu;
  double nonmatch_initial_value = -penalties->log_H + penalties->mu;

  double active_value = cohort->match_is_rare(site_index, a) ? match_initial_value : nonmatch_initial_value;
  double default_value = cohort->match_is_rare(site_index, a) ? nonmatch_initial_value : match_initial_value;
  
  std::fill(R.begin(), R.end(), default_value);
  
  if(cohort->number_active(site_index, a) != 0) {
    const rowSet& active_rows = cohort->get_active_rowSet(site_index, a);
    rowSet::const_iterator it = active_rows.begin();
    rowSet::const_iterator rows_end = active_rows.end();
    for(it; it != rows_end; ++it) {
      R[*it] = active_value;
    }
  }

  if(cohort->number_matching(site_index, a) == 0) {
    S = penalties->mu;
  } else if(cohort->number_not_matching(site_index, a) == 0) {
    S = penalties->one_minus_mu;
  } else {  
    S = -penalties->log_H + 
                logsum(log(cohort->number_matching(site_index, a)) + penalties->one_minus_mu,
                       log(cohort->number_not_matching(site_index, a)) + penalties->mu);
  }
  record_last_extended(a);
}

void fastFwdAlgState::update_subset_of_Rs(const rowSet& indices,
              bool active_is_match) {
  double correction = penalties->get_minority_map_correction(active_is_match);
  rowSet::const_iterator it = indices.begin();
  rowSet::const_iterator rows_end = indices.end();
  for(it; it != rows_end; ++it) {
    size_t row = *it;
    R[row] = correction + calculate_R(R[row], map.get_map(row));
  }
}

void fastFwdAlgState::fast_update_S(const rowSet& indices, bool active_is_match) {
  penalties->update_S(S, R, indices.begin(), indices.end(), active_is_match);
}

void fastFwdAlgState::extend_probability_at_site(size_t site_index,
            alleleValue a) {
  bool match_is_rare = cohort->match_is_rare(site_index, a);
  DPUpdateMap current_map = penalties->get_current_map(S, match_is_rare);
  const rowSet& active_rows = cohort->get_active_rowSet(site_index, a);
  extend_probability_at_site(current_map, active_rows, match_is_rare, a);
}

void fastFwdAlgState::extend_probability_at_span_after(size_t site_index,
            size_t mismatch_count = 0) {
  size_t length = reference->span_length_after(site_index);
  extend_probability_at_span_after_anonymous(length, mismatch_count);
}

void fastFwdAlgState::take_snapshot() {
  if(last_extended >= 0) {
    map.hard_update_all();
    size_t j = last_extended;
    bool reference_is_homogeneous = (cohort->number_matching(j, last_allele) == 0 || cohort->number_not_matching(j, last_allele) == 0);
    if(reference_is_homogeneous || last_extended_is_span()) {
      for(size_t i = 0; i < cohort->get_n_haplotypes(); i++) {
        R[i] = calculate_R(R[i], map.get_map(i));
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
          R[i] = calculate_R(R[i], map.get_map(i));
        }
      }
    }
    // since all R-values are up to date, we do not need entries in the lazyEvalMap
    // therefore we can clear them all and replace them with the identity map
    map.hard_clear_all();
  }
} 

double fastFwdAlgState::prefix_likelihood() const {
  return S;
}

double fastFwdAlgState::partial_likelihood_by_row(size_t row) const {
  return R[row];
}

double calculate_R(double oldR, const DPUpdateMap& map) {
  return map.of(oldR);
}

double calculate_R(double oldR, double coefficient, double constant) {
  return calculate_R(oldR, DPUpdateMap(coefficient, constant));
}

void fastFwdAlgState::extend_probability_at_site(const DPUpdateMap& current_map, 
            const rowSet& active_rows, bool match_is_rare, 
            alleleValue a) {
  map.stage_map_for_site(current_map);
  if(active_rows.empty() && match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->mu + S;
  } else if(active_rows.empty() && !match_is_rare) {
    // separate case to avoid log-summing "log 0"
    S = penalties->one_minus_mu + S;
  } else {
    map.update_active_rows(active_rows);
    update_subset_of_Rs(active_rows, match_is_rare);
    fast_update_S(active_rows, match_is_rare);
    map.reset_rows(active_rows);
  }
  record_last_extended(a);
  return;
}

void fastFwdAlgState::extend_probability_at_site(
            const rowSet& active_rows, bool match_is_rare, 
            alleleValue a) {
  DPUpdateMap current_map = penalties->get_current_map(S, match_is_rare);
  extend_probability_at_site(current_map, active_rows, match_is_rare, a);
}

void fastFwdAlgState::extend_probability_at_span_after_anonymous(size_t l, size_t mismatch_count) {
  double m = penalties->span_mutation_penalty(l, mismatch_count);
  double coefficient = penalties->pow_rho_c(l);
  double constant = penalties->span_polynomial(l) + S;
  map.stage_map_for_span(DPUpdateMap(m + coefficient, constant - coefficient));
  S = m + S;
  last_span_extended = last_extended;
}

slowFwdSolver::slowFwdSolver(siteIndex* ref, const penaltySet* pen, const haplotypeCohort* haplotypes) :
            reference(ref), penalties(pen), cohort(haplotypes) {
}
            
double slowFwdSolver::calculate_probability_quadratic(const vector<alleleValue>& q, size_t start_site = 0) {
  R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H);
  for(size_t j = 0; j < R.size(); j++) {
    bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    R[j] += emission;
  }
  vector<double> last_R;
  double same_transition = log1p(-exp(penalties->rho) * (penalties->H - 1));
  for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
    last_R = R;
    if(reference->has_span_after(i - 1)) {
      size_t l = reference->span_length_after(i - 1);
      double coefficient = penalties->span_mutation_penalty(l, 0) + penalties->composed_R_coefficient(l);
      double constant = penalties->span_coefficient(l) + S - penalties->composed_R_coefficient(l);
      for(size_t j = 0; j < R.size(); j++) {
        R[j] = coefficient + logsum(R[j], constant);
      }
    }
    S = log_big_sum(R);
    for(size_t j = 0; j < R.size(); j++) {
      vector<double> temp = last_R;
      for(size_t k = 0; k < temp.size(); k++) {
        temp[k] += j == k ? same_transition : penalties->rho;
      }
      R[j] = log_big_sum(temp);
      bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      R[j] += emission;
    }
  }
  S = log_big_sum(R);
  return S;
}

double slowFwdSolver::calculate_probability_linear(const vector<alleleValue>& q, size_t start_site = 0) {
  R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H);
  for(size_t j = 0; j < R.size(); j++) {
    bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    R[j] += emission;
  }
  vector<double> last_R;
  for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
    last_R = R;
    S = log_big_sum(last_R);
    // do span stuff
    if(reference->has_span_after(i - 1)) {
      size_t l = reference->span_length_after(i - 1);
      double coefficient = penalties->span_mutation_penalty(l, 0) + penalties->composed_R_coefficient(l);
      double constant = penalties->span_coefficient(l) + S - penalties->composed_R_coefficient(l);
      for(size_t j = 0; j < R.size(); j++) {
        R[j] = coefficient + logsum(R[j], constant);
      }
    }
    S = log_big_sum(R);
    for(size_t j = 0; j < R.size(); j++) {
      R[j] = logsum(S + penalties->rho, last_R[j] + penalties->R_coefficient); 
      bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      R[j] += emission;
    }
  }
  S = log_big_sum(R);
  return S;
}

double slowFwdSolver::calculate_probability_quadratic(const inputHaplotype* q) {
  vector<alleleValue> alleles = q->get_alleles();
  size_t start_site = q->get_start_site();
  S = q->has_left_tail() ? penalties->span_mutation_penalty(q->get_left_tail(), 0) : 0;
  R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H + S);
  for(size_t j = 0; j < R.size(); j++) {
    bool matches = (cohort->allele_at(start_site, j) == alleles[0]);
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    R[j] += emission;
  }
  S = log_big_sum(R);

  vector<double> last_R;
  double same_transition = log1p(-exp(penalties->rho) * (penalties->H - 1));
  for(size_t i = 0; i < q->number_of_sites(); i++) {
    // do span stuff
    if(q->has_span_after(i - 1)) {
      size_t l = q->get_span_after(i - 1);
      double coefficient = penalties->pow_rho_c(l);
      double constant = penalties->span_polynomial(l) + S;
      for(size_t j = 0; j < R.size(); j++) {
        R[j] = penalties->span_mutation_penalty(l, 0) + logsum(coefficient + R[j], constant);
      }
    }
    S = log_big_sum(R);
    last_R = R;
    for(size_t j = 0; j < R.size(); j++) {
      vector<double> temp = last_R;
      for(size_t k = 0; k < temp.size(); k++) {
        temp[k] += j == k ? same_transition : penalties->rho;
      }
      R[j] = log_big_sum(temp);
      bool matches = (cohort->allele_at(start_site + i, j) == alleles[i]);
      double emission = matches ? penalties->one_minus_mu : penalties->mu;
      R[j] += emission;
    }
    S = log_big_sum(R);
  }
  if(q->has_span_after(q->number_of_sites() - 1)) {
    size_t l = q->get_span_after(q->number_of_sites() - 1);
    double coefficient = penalties->pow_rho_c(l);
    double constant = penalties->span_polynomial(l) + S;
    for(size_t j = 0; j < R.size(); j++) {
      R[j] = penalties->span_mutation_penalty(l, 0) + logsum(coefficient + R[j], constant);
    }
  }
  
  S = log_big_sum(R);
  return S;
}

// double slowFwdSolver::calculate_probability_linear(const inputHaplotype* q) {
//   vector<alleleValue> alleles = q->get_alleles();
//   size_t start_site = q->get_start_site();
//   S = q->has_left_tail() ? penalties->span_mutation_penalty(q->get_left_tail(), 0) : 0;
//   R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H + S);
//   for(size_t j = 0; j < R.size(); j++) {
//     bool matches = (cohort->allele_at(start_site, j) == alleles[0]);
//     double emission = matches ? penalties->one_minus_mu : penalties->mu;
//     R[j] += emission;
//   }
//   S = log_big_sum(R);
// 
//   for(size_t i = 0; i < q->number_of_sites(); i++) {
//     // do span stuff
//     if(q->has_span_after(i - 1)) {
//       size_t l = q->get_span_after(i - 1);
//       double coefficient = penalties->pow_rho_c(l);
//       double constant = penalties->span_polynomial(l) + S - penalties->log_H;
//       for(size_t j = 0; j < R.size(); j++) {
//         R[j] = penalties->span_mutation_penalty(l, 0) + logsum(coefficient + R[j], constant);
//       }
//     }
//     S = log_big_sum(R);
//     for(size_t j = 0; j < R.size(); j++) {
//       R[j] = logsum(S + penalties->rho, R[j] + penalties->R_coefficient); 
//       bool matches = (cohort->allele_at(start_site + i, j) == alleles[i]);
//       double emission = matches ? penalties->one_minus_mu : penalties->mu;
//       R[j] += emission;
//     }
//     S = log_big_sum(R);
//   }
//   if(q->has_span_after(q->number_of_sites() - 1)) {
//     size_t l = q->get_span_after(q->number_of_sites() - 1);
//     double coefficient = penalties->pow_rho_c(l);
//     double constant = penalties->span_polynomial(l) + S - penalties->log_H;
//     for(size_t j = 0; j < R.size(); j++) {
//       R[j] = penalties->span_mutation_penalty(l, 0) + logsum(coefficient + R[j], constant);
//     }
//   }  
//   
//   S = log_big_sum(R);
//   return S;
// }

double slowFwdSolver::calculate_probability_linear(const inputHaplotype* q) {
  initialize_linear(q);
  for(size_t i = 0; i < q->number_of_sites(); i++) {
    if(q->has_span_after(i - 1)) {
      extend_span_linear(q, i - 1);
    }
    extend_site_linear(q, i);
  }
  if(q->has_span_after(q->number_of_sites() - 1)) {
    extend_span_linear(q, q->number_of_sites() - 1);
  }
  return S;
}


void slowFwdSolver::initialize_linear(const inputHaplotype* q) {
  size_t start_site = q->get_start_site();
  S = q->has_left_tail() ? penalties->span_mutation_penalty(q->get_left_tail(), 0) : 0;
  R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H + S);
  for(size_t j = 0; j < R.size(); j++) {
    bool matches = (cohort->allele_at(start_site, j) == q->get_allele(0));
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    R[j] += emission;
  }
  S = log_big_sum(R);
}

void slowFwdSolver::extend_site_linear(const inputHaplotype* q, size_t site) {
  size_t start_site = q->get_start_site();
  for(size_t j = 0; j < R.size(); j++) {
    R[j] = logsum(S + penalties->rho, R[j] + penalties->R_coefficient); 
    bool matches = (cohort->allele_at(start_site + site, j) == q->get_allele(site));
    double emission = matches ? penalties->one_minus_mu : penalties->mu;
    R[j] += emission;
  }
  // S = log_big_sum(R);
  if(!cohort->match_is_rare(start_site + site, q->get_allele(site))) {
    if(cohort->number_not_matching(start_site + site, q->get_allele(site)) != 0) {
      const rowSet& indices = cohort->get_active_rowSet(start_site + site, q->get_allele(site));
      rowSet::const_iterator it = indices.begin();
      rowSet::const_iterator rows_end = indices.end();
      // vector<double> to_correct;
      // for(it; it != rows_end; ++it) {
      //   size_t row = *it;
      //   to_correct.push_back(R[row]);
      // }
      S = logdiff(penalties->one_minus_mu + S, penalties->one_minus_2mu - penalties->mu + log_big_sum(it, rows_end, R));
    } else {
      S += penalties->one_minus_mu;
    }
  } else {
    if(cohort->number_matching(start_site + site, q->get_allele(site)) != 0) {
      const rowSet& indices = cohort->get_active_rowSet(start_site + site, q->get_allele(site));
      rowSet::const_iterator it = indices.begin();
      rowSet::const_iterator rows_end = indices.end();
      // vector<double> to_correct;
      // for(it; it != rows_end; ++it) {
      //   size_t row = *it;
      //   to_correct.push_back(R[row]);
      // }
      S = logsum(penalties->mu + S, penalties->one_minus_2mu - penalties->one_minus_mu + log_big_sum(it, rows_end, R));
    } else {
      S += penalties->mu;
    }
  }
}

// void slowFwdSolver::extend_site_linear(const inputHaplotype* q, size_t site) {
//   size_t start_site = q->get_start_site();
//   for(size_t j = 0; j < R.size(); j++) {
//     R[j] = logsum(S + penalties->rho, R[j] + penalties->R_coefficient); 
//     bool matches = (cohort->allele_at(start_site + site, j) == q->get_allele(site));
//     double emission = matches ? penalties->one_minus_mu : penalties->mu;
//     R[j] += emission;
//   }
//   S = log_big_sum(R);
// }

void slowFwdSolver::extend_span_linear(const inputHaplotype* q, size_t site) {
  if(q->has_span_after(site)) {
    size_t l = q->get_span_after(site);
    double coefficient = penalties->pow_rho_c(l);
    double constant = penalties->span_polynomial(l) + S;
    for(size_t j = 0; j < R.size(); j++) {
      R[j] = penalties->span_mutation_penalty(l, 0) + logsum(coefficient + R[j], constant);
    }
    S += penalties->span_mutation_penalty(l, 0);
    // S = log_big_sum(R);
  }
}

// pair<vector<double>, vector<size_t> > slowFwdSolver::sequence_statistics(const vector<alleleValue>& q, size_t start_site) {
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
//   vector<size_t> active_rows(q.size(), 0);
//   vector<size_t> history_length(q.size(), 0);
//   vector<size_t> active_history_length(q.size(), 0);
//   vector<size_t> map_classes(q.size(), 0);
//   vector<size_t> active_map_classes(q.size(), 0);
//   
//   vector<size_t> times_active(cohort->get_n_haplotypes(), 
//   0);
//   vector<size_t> last_updated(cohort->get_n_haplotypes(), start_site);
//   R = vector<double>(cohort->get_n_haplotypes(), -penalties->log_H);
//   for(size_t j = 0; j < R.size(); j++) {
//     bool matches = (cohort->allele_at(0, j) == q[0 - start_site]);
//     double emission = matches ? penalties->one_minus_mu : penalties->mu;
//     R[j] += emission;
//   }
//   vector<double> last_R;
//   for(size_t i = start_site + 1; i < start_site + q.size(); i++) {
//     last_R = R;
//     S = log_big_sum(last_R);
//     // do span stuff
//     if(reference->has_span_after(i - 1)) {
//       size_t l = reference->span_length_after(i - 1);
//       double coefficient = penalties->span_mutation_penalty(l, 0) + penalties->composed_R_coefficient(l);
//       double constant = penalties->span_coefficient(l) - penalties->composed_R_coefficient(l);
//       for(size_t j = 0; j < R.size(); j++) {
//         R[j] = coefficient + logsum(R[j], constant);
//       }
//     }
//     for(size_t j = 0; j < R.size(); j++) {
//       R[j] = logsum(S + penalties->rho, last_R[j] + penalties->R_coefficient); 
//       bool matches = (cohort->allele_at(i, j) == q[i - start_site]);
//       double emission = matches ? penalties->one_minus_mu : penalties->mu;
//       R[j] += emission;
//     }
//     S = log_big_sum(R);
//     // DATA gathering  
// 
//     double eps = 0.000000001;
//     // VALUES --------------------------------------------------------------
//     for(size_t j = 0; j < cohort->get_n_haplotypes(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(fabs(R.at(j) - R.at(k)) < eps) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++unique_values.at(i - start_site);
//       }
//     }
//     vector<size_t> i_active_rows = cohort->get_active_rows(i, q.at(i - start_site));
//     for(size_t j = 0; j < i_active_rows.size(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(fabs(R.at(i_active_rows.at(j)) - R.at(i_active_rows.at(k))) < eps) {
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
//     active_rows.at(i - start_site) = i_active_rows.size();
//     // HISTORY LENGTH -----------------------------------------------------------
//     size_t oldest_active = i;
//     size_t n_active_mapclasses = 0;
//     size_t oldest = i;
//     size_t n_mapclasses = 0;
//     for(size_t j = 0; j < i_active_rows.size(); j++) {
//       bool found = false;
//       for(size_t k = 0; k < j; k++) {
//         if(last_updated.at(i_active_rows.at(j)) == last_updated.at(i_active_rows.at(k))) {
//           found = true;
//           break;
//         }
//       }
//       if(!found) {
//         ++n_active_mapclasses;
//       }
//       if(last_updated.at(i_active_rows.at(j)) < oldest_active) {
//         oldest_active = last_updated.at(i_active_rows.at(j));
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
//     for(size_t j = 0; j < i_active_rows.size(); j++) {
//       last_updated.at(i_active_rows.at(j)) = i;
//       ++times_active.at(i_active_rows.at(j));
//     }
//   }
//   
//   size_t max_u_value = 1;
//   double avg_u_value = 0;
//   size_t max_u_active_value = 1;
//   double avg_u_active_value = 0;
//   size_t max_active_rows = 0;
//   double avg_active_rows = 0;
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
//     if(active_rows.at(i) > max_active_rows) {
//       max_active_rows = active_rows.at(i);
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
//     avg_active_rows += active_rows.at(i);
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
//   avg_active_rows = avg_active_rows/common_denom;
//   avg_history_length = avg_history_length/common_denom;
//   avg_active_history_length = avg_active_history_length/common_denom;
//   avg_map_classes = avg_map_classes/common_denom;
//   avg_active_map_classes = avg_active_map_classes/common_denom;
//   avg_times_active = avg_times_active/times_active.size();
//   
//   vector<size_t> to_return_max = {
//     max_u_value,
//     max_u_active_value,
//     max_active_rows,
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
//     avg_active_rows,
//     avg_history_length,
//     avg_active_history_length,
//     avg_map_classes,
//     avg_active_map_classes,
//     avg_times_active
//   };
//   
//   return make_pair(to_return_avg, to_return_max);
// }
