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
  double ft_of_one;
  // log of (1 - 2*rho + H*rho)
  double fs_of_one;
  
  penaltySet(double logRho, double logMu, int H);
  ~penaltySet();
};

// A haplotypeMatrix is the matrix which iteratively calculates haplotype
// likelihood. It takes in a haplotypeCohort and a linearReferenceStructure, an
// inputHaplotype built against the linearReferenceStructure, and a penaltySet
// of mutation and recombination penalties. It calculates and returns the
// likelihood of the inputHaplotype relative to the haplotypeCohort when
// calculate_probability is called
struct haplotypeMatrix{
private:
  linearReferenceStructure* reference;
  haplotypeCohort* cohort;
  penaltySet* penalties;
  delayMap map;
  
  // trackers for the last indices extended. spans are indexed according to
  // the site preceding them, except for index -1, the span before site 0
  // -1 : nothing extended; i : index i last extended
  int last_extended = -1;
  size_t last_site_extended;
  alleleValue last_allele;
  void record_last_extended(alleleValue a);
  // -2 : nothing extended; indexing as above
  int last_span_extended = -2;
  bool last_extended_is_span();

public:
  haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
            haplotypeCohort* haplotypes);
  haplotypeMatrix(const haplotypeMatrix& other, bool copy_map);
  ~haplotypeMatrix();
  
  double S;
  vector<double> R;
  
  delayMap get_map();

  // TODO: implement these 
  // double estimate_probability_at_site(size_t j, alleleValue a);
  // helper functions for probability estimation
  // double minContinuing(int j);
  // double minMutating(int j);
  // double maxSwitching(int j);
  
  // checks for a span before the first site
  void initialize_probability(inputHaplotype* q);
  void extend_probability_at_site(inputHaplotype* q, size_t j);
  void extend_probability_at_span_after(inputHaplotype* q, size_t j);
  
  double calculate_probability(inputHaplotype* q);
  
  void initialize_probability(size_t site_index, alleleValue a,
              size_t left_tail_length = 0, size_t augmentation_count = 0);
  void initialize_probability_at_span(size_t length, size_t augmentation_count);
  void initialize_probability_at_site(size_t site_index, alleleValue a);
  void extend_probability_at_site(size_t site_index, alleleValue a);
  void extend_probability_at_span_after(size_t site_index, 
              size_t augmentation_count);            
              
  void take_snapshot();
  double prefix_likelihood();
  double partial_likelihood_by_row(size_t row);
};

double calculate_R(double oldR, DPUpdateMap map);
double calculate_R(double oldR, double coefficient, double constant);

#endif