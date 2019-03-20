#ifndef LINEAR_HAPLO_PROBABILITY_H
#define LINEAR_HAPLO_PROBABILITY_H

#include "reference.hpp"
#include "input_haplotype.hpp"
#include "DP_map.hpp"
#include "delay_multiplier.hpp"
#include "penalty_set.hpp"

using namespace std;

struct liStephensModel{
  liStephensModel(siteIndex* reference, haplotypeCohort* cohort, const penaltySet* penalties);
  siteIndex* reference;
  haplotypeCohort* cohort;
  const penaltySet* penalties;
};

// A fastFwdAlgState contains the DP matrix for our forward algorithm analogue.
// It takes in a haplotypeCohort and a siteIndex, an
// inputHaplotype built against the siteIndex, and a penaltySet
// of mutation and recombination penalties. It calculates and returns the
// likelihood of the inputHaplotype relative to the haplotypeCohort when
// calculate_probability is called
struct fastFwdAlgState{
private:
//-- underlying structures -----------------------------------------------------
  
  const siteIndex* reference;
  const haplotypeCohort* cohort;
  const penaltySet* penalties;
  
//-- lazy evaluation "backend" -------------------------------------------------
  
  delayedEvalMap map;
  
//-- position markers ----------------------------------------------------------

  // trackers for the last indices extended
  typedef int64_t dp_column_t;
  const static dp_column_t SITE_UNEXTENDED;
  const static dp_column_t INITIAL_SPAN;
  const static dp_column_t SPAN_UNEXTENDED;
  
  dp_column_t last_extended;
  dp_column_t last_span_extended;
  alleleValue last_allele;
  
  double sum;
  vector<double> rows;
  
  inline void record_last_extended(alleleValue a);

#ifdef DEBUG
public:
#endif
  inline bool last_extended_is_span() const;
  inline size_t get_last_site() const;
  
  void initialize_probability(const inputHaplotype* observed);
  void initialize_probability_at_span(size_t length);
  void initialize_probability_at_site(size_t site_index, alleleValue a);
  
  void update_sum(const rowSet& indices, bool active_is_match);
  
  void extend_probability_at_site(const DPUpdateMap& current_map, 
                                  const rowSet& active_row_indices, 
                                  bool match_is_rare, 
                                  alleleValue a);
  void extend_probability_at_site(size_t i, alleleValue a);
  void extend_probability_at_site(alleleValue a);
  
  void extend_probability_at_span_after_abstract(size_t l);
  void extend_probability_at_span();
  
#ifdef DEBUG  
  double get_sum() const;
  double get_row() const;
#endif
  void take_snapshot();
  
#ifndef DEBUG
public:
#endif
  fastFwdAlgState(siteIndex* ref, const penaltySet* pen,
    const haplotypeCohort* haplotypes);
  fastFwdAlgState(const fastFwdAlgState& other, bool copy_map);
  ~fastFwdAlgState();

  double calculate_probability(const inputHaplotype* observed);
  
#ifdef TIME_PROBABILITY_INTERNALS
  double* t_total;
  double* t_readwrite;
  double* t_delay;
  void set_timers(double* t_total, double* t_readwrite, double* t_delay);
#endif
};

struct slowFwdSolver{
  siteIndex* reference;
  const penaltySet* penalties;
  const haplotypeCohort* cohort;
  vector<double> rows;
  double sum;
  slowFwdSolver(siteIndex* ref, const penaltySet* pen, const haplotypeCohort* haplotypes);
  double calculate_probability_quadratic(const inputHaplotype* observed_haplotype);
  double calculate_probability_linear(const inputHaplotype* observed_haplotype);
  
  void initialize_linear(const inputHaplotype* q);
  void extend_site_linear(const inputHaplotype* q, size_t site);
  void extend_span_linear(const inputHaplotype* q, size_t site);
};

double calculate_row(double old_row, const DPUpdateMap& map);
double calculate_row(double old_row, double coefficient, double constant);

#endif