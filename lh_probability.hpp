#ifndef LINEAR_HAPLO_PROBABILITY_H
#define LINEAR_HAPLO_PROBABILITY_H

#include "lh_math.hpp"
#include "lh_reference.hpp"
#include "lh_input_haplotype.hpp"
#include "lh_DP_map.hpp"
#include "lh_delay_multiplier.hpp"

using namespace std;

// stores a shared set of penalty-derived coefficients for use in calculations
// according to our model
struct penaltySet{
  int H;
  double log_H;
  double rho;
  double mu;
  double one_minus_rho;
  double one_minus_mu;
  double one_minus_2mu;
  
  // log of (1 - 2*rho)
  double alpha_value;
  // log of (1 - 2*rho + H*rho)
  double beta_value;
  
  penaltySet(double logRho, double logMu, int H);
  ~penaltySet();
  
  double span_coefficient(size_t l) const;
  double alpha(size_t l) const;
  double beta(size_t l) const;
  double span_mutation_penalty(size_t l, size_t a) const;
  
  DPUpdateMap get_match_map(double last_sum) const;
  DPUpdateMap get_non_match_map(double last_sum) const;
  DPUpdateMap get_current_map(double last_sum, bool match_is_rare) const;
  double get_minority_map_correction(bool match_is_rare) const;
  void update_S(double& S, const vector<double>& summands, bool match_is_rare) const;
};

// A haplotypeMatrix is the matrix which iteratively calculates haplotype
// likelihood. It takes in a haplotypeCohort and a linearReferenceStructure, an
// inputHaplotype built against the linearReferenceStructure, and a penaltySet
// of mutation and recombination penalties. It calculates and returns the
// likelihood of the inputHaplotype relative to the haplotypeCohort when
// calculate_probability is called
struct haplotypeMatrix{
private:
  const linearReferenceStructure* reference;
  const haplotypeCohort* cohort;
  const penaltySet* penalties;
  delayMap map;
  
  // trackers for the last indices extended. spans are indexed according to
  // the site preceding them, except for index -1, the span before site 0
  // -1 : nothing extended; i : index i last extended
  int last_extended = -1;
  // -2 : nothing extended; indexing as above
  int last_span_extended = -2;
  
  alleleValue last_allele;
  void record_last_extended(alleleValue a);
  
public:
  haplotypeMatrix(const linearReferenceStructure* ref, const penaltySet* pen,
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
};

double calculate_R(double oldR, const DPUpdateMap& map);
double calculate_R(double oldR, double coefficient, double constant);

#endif