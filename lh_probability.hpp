#ifndef LINEAR_HAPLO_PROBABILITY_H
#define LINEAR_HAPLO_PROBABILITY_H

#include "lh_reference.hpp"
#include "lh_input_haplotype.hpp"

using namespace std;

// stores a shared set of penalty-derived coefficients for use in calculations
// according to our model
struct penaltySet{
  int H;
  double log_H;
  double log_rho;
  double log_mu;
  double log_rho_complement;
  double log_mu_complement;
  double log_2mu_complement;
  
  // log of (1 - 2*rho)
  double log_ft_base;
  // log of (1 - 2*rho + H*rho)
  double log_fs_base;
  
  penaltySet(double logRho, double logMu, int H);
  ~penaltySet();
};

struct haplotypeMatrix{
private:
  linearReferenceStructure* reference;
  haplotypeCohort* cohort;
  penaltySet* penalties;
  inputHaplotype* query;
  
  // This is used when the inputHaplotype begins with a span
  double initial_R;
  
  // trackers for the last indices extended. spans are indexed according to
  // the site preceding them, except for index -1, the span before site 0
  // -1 : nothing extended; i : index i last extended
  int last_extended = -1;
  // -2 : nothing extended; indexing as above
  int last_span_extended = -2;
  
  bool last_extended_is_span();
  
  double last_S();
  double last_R(size_t index_among_haplotypes);
  
  size_t size();
  
  void initialize_probability();
  void extend_probability_at_site(size_t j, alleleValue a);
  void extend_probability_at_site(size_t j);
  void extend_probability_at_span_after(size_t j, 
              int augmentation_count);
public:
  vector<double> S;
  vector<vector<double> > R;
  
  haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
            haplotypeCohort* haplotypes, inputHaplotype* query);
  ~haplotypeMatrix();
  double calculate_probabilities();

  // TODO: implement these 
  // double estimate_probability_at_site(size_t j, alleleValue a);
  // helper functions for probability estimation
  // double minContinuing(int j);
  // double minMutating(int j);
  // double maxSwitching(int j);
  
  // vector<size_t> get_matches(size_t index, alleleValue a);
  // size_t number_matching(size_t index, alleleValue a);
  // size_t number_not_matching(size_t index, alleleValue a);
};

// log-space math functions
double logdiff(double a, double b);
double logsum(double a, double b);
double log_big_sum(vector<double>& logRs);
double log_weighted_big_sum(vector<double>& logRs, vector<int>& counts);

#endif