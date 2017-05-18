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
  // log of ((1 - rho) + (H - 1)*rho)
  double log_fs_base;
  // log of (1-mu)(1-2*rho)
  double log_R_match_coefficient;
  double log_R_mismatch_coefficient;
  // log of (1-mu)*rho
  double log_match_rho;
  double log_mismatch_rho;
  
  penaltySet(double logRho, double logMu, int H);
  ~penaltySet();
};

struct haplotypeMatrix{
private:
  linearReferenceStructure* reference;
  haplotypeCohort* cohort;
  penaltySet* penalties;
  inputHaplotype* query;
    
  vector<double> initial_R;
  
  // trackers for the last indices extended. spans are indexed according to
  // the site preceding them, except for index -1, the span before site 0
  // -1 : nothing extended; i : index i last extended
  int last_extended = -1;
  // -2 : nothing extended; indexing as above
  int last_span_extended = -2;
  
  bool last_extended_is_span();
  
  double last_S();
  double last_R(int i);
  
  size_t size();
  
  void initialize_probability();
  
  void extend_probability_at_site(int j, alleleValue a);
  void extend_probability_at_site(int j);
  void extend_probability_at_span_after(int j, int augmentation_count);
public:
  vector<double> S;
  vector<vector<double> > R;
  
  haplotypeMatrix(linearReferenceStructure* ref, penaltySet* pen,
            haplotypeCohort* haplotypes, inputHaplotype* query);
  ~haplotypeMatrix();
  double calculate_probabilities();

  // TODO: implement these 
  // double estimate_probability_at_site(int j, alleleValue a);
  // helper functions for probability estimation
  // double minContinuing(int j);
  // double minMutating(int j);
  // double maxSwitching(int j);
};

// log-space math functions
double logdiff(double a, double b);
double logsum(double a, double b);
double log_big_sum(vector<double>& logRs);
double log_weighted_big_sum(vector<double>& logRs, vector<int>& counts);

#endif