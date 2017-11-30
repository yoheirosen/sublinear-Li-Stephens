#ifndef LINEAR_HAPLO_PROBABILITY_H
#define LINEAR_HAPLO_PROBABILITY_H

#include "math.hpp"
#include "reference.hpp"
#include "input_haplotype.hpp"
#include "DP_map.hpp"
#include "delay_multiplier.hpp"
#include "penalty_set.hpp"

using namespace std;

// A haplotypeMatrix is the matrix which iteratively calculates haplotype
// likelihood. It takes in a haplotypeCohort and a siteIndex, an
// inputHaplotype built against the siteIndex, and a penaltySet
// of mutation and recombination penalties. It calculates and returns the
// likelihood of the inputHaplotype relative to the haplotypeCohort when
// calculate_probability is called
struct haplotypeMatrix{
private:
  const siteIndex* reference;
  const haplotypeCohort* cohort;
  const penaltySet* penalties;
  // TODO (10/3/17) streamline memory use here
  delayMap map;
  
    // trackers for the last indices extended
  // -1 : nothing extended; i : index i last extended
  int last_extended = -1;
  // -2 : nothing extended; -1 : before-first-site span extended; i : span after
  // index i extended
  int last_span_extended = -2;
  alleleValue last_allele;
    void record_last_extended(alleleValue a);
  
public:
  haplotypeMatrix(const siteIndex* ref, const penaltySet* pen,
            const haplotypeCohort* haplotypes);
  haplotypeMatrix(const haplotypeMatrix& other, bool copy_map);
  ~haplotypeMatrix();
  
  double S;
  vector<double> R;
  
  delayMap& get_maps();
  
  // checks for a span before the first site
  void initialize_probability(const inputHaplotype* q);
  void extend_probability_at_site(const inputHaplotype* q, size_t j);
  void extend_probability_at_span_after(const inputHaplotype* q, size_t j);
  
  void extend_probability_at_site(const DPUpdateMap& current_map, 
              const vector<size_t>& active_rows, bool match_is_rare, 
              alleleValue a);
  void extend_probability_at_site(const vector<size_t>& active_rows, 
              bool match_is_rare, alleleValue a);
  void extend_probability_at_span_after_anonymous(size_t l,
              size_t mismatch_count);
  
  double calculate_probability(const inputHaplotype* q);
  
  void initialize_probability(size_t site_index, alleleValue a,
              size_t left_tail_length = 0, size_t mismatch_count = 0);
  void initialize_probability_at_span(size_t length, size_t mismatch_count);
  void initialize_probability_at_site(size_t site_index, alleleValue a);
  void extend_probability_at_site(size_t site_index, alleleValue a);
  void extend_probability_at_span_after(size_t site_index, 
              size_t mismatch_count);            
  
  void update_subset_of_Rs(const vector<size_t>& indices, bool active_is_match);
  void fast_update_S(const vector<size_t>& indices, bool active_is_match);
  
  // Applies the delayed-arithmetic maps to all R-values 
  // This also calls hard_clear_all on the delayed-arithmetic map, which
  // resets all linear maps contained to the identity map; this means that
  // accidentally calling take_snapshot() twice has no effect
  void take_snapshot();
  double prefix_likelihood() const;
  double partial_likelihood_by_row(size_t row) const;
  
  bool last_extended_is_span() const;
  size_t get_last_site() const;
  
  void extend_probability_at_site(const DPUpdateMap& current_map, 
              const rowSet& active_rows, bool match_is_rare, 
              alleleValue a);
  void extend_probability_at_site(const rowSet& active_rows, 
              bool match_is_rare, alleleValue a);
  void update_subset_of_Rs(const rowSet& indices, bool active_is_match);
  void fast_update_S(const rowSet& indices, bool active_is_match);
};

struct forward_solver{
  const siteIndex* reference;
  const penaltySet* penalties;
  const haplotypeCohort* cohort;
  forward_solver(const siteIndex* ref, const penaltySet* pen,
            const haplotypeCohort* haplotypes);
  double calculate_probability_quadratic(const vector<alleleValue>& q, size_t start_site);
  double calculate_probability_linear(const vector<alleleValue>& q, size_t start_site);
};

double calculate_R(double oldR, const DPUpdateMap& map);
double calculate_R(double oldR, double coefficient, double constant);

#endif