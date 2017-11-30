// A population reference cohort is expressed using two elements:
// 1. a linearReferenceStructure which stores positions of sites of variation
//    in a population
// 2. a haplotypeCohort which stores vectors of haplotype indices containing a
//    given allele at a given site. This is essential for rapid calculation of
//    the population likelihoods

#ifndef LINEAR_REFERENCE_STRUCTURE_H
#define LINEAR_REFERENCE_STRUCTURE_H

#include <string>
#include <vector>
#include <unordered_map>
#include "allele.hpp"
#include "row_set.hpp"

using namespace std;

struct linearReferenceStructure{
private:
  size_t global_offset = 0;
  size_t length = 0;
  bool final_span_calculated = false;
  bool alleles_set = false;
  unordered_map<size_t, size_t> position_to_site_index;
  vector<size_t> site_index_to_position;
  vector<alleleValue> site_index_to_reference_allele;
  
  void add_site(size_t position, alleleValue reference_value);
  
  vector<size_t> span_lengths;
  size_t leading_span_length;
public:
  linearReferenceStructure(const vector<string>& haplotypes,
            const string& reference_values);
  linearReferenceStructure(const vector<size_t>& positions, size_t length,
            const vector<alleleValue>& reference_values);
  linearReferenceStructure(const vector<size_t>& positions, size_t length, size_t global_offset);
  linearReferenceStructure(const vector<size_t>& positions,
            const string& reference_sequence);
  linearReferenceStructure(size_t global_offset);
  ~linearReferenceStructure();
  
  bool is_site(size_t actual_position) const;
  size_t get_site_index(size_t actual_position) const;
  size_t get_position(size_t site_index) const;
  
  bool has_span_before(size_t site_index) const;
  bool has_span_after(size_t site_index) const;
  size_t span_length_before(size_t site_index) const;
  size_t span_length_after(size_t site_index) const;
  
  size_t number_of_sites() const;
  size_t absolute_length();
  size_t absolute_length() const;
  
  alleleValue get_reference_allele_at_site(size_t site_index) const;
  
  // behavior when there is no site above: returns the 1-past-the-end index
  size_t find_site_above(size_t position) const;
  // behavior when there is no site below: returns SIZE_MAX
  size_t find_site_below(size_t position) const;
  
  // >= 0 : successful, returns site-index
  // -1 : out of order
  // -2 : collision
  // -3 : site-vector locked
  int64_t add_site(size_t position);
  void set_initial_span(size_t length);
  
  void calculate_final_span_length(size_t reference_length);
  
  size_t pos_ref2global(size_t p) const;
  size_t start_position() const;
  int64_t pos_global2ref(int64_t p) const;
  
  // 1: success
  // 0: not a site
  // -1: out of bounds
  int set_allele_at_global_pos(size_t p, alleleValue allele);
  bool set_allele_at_pos(size_t p, alleleValue allele);
  void set_allele_at_site(size_t site, alleleValue allele);
};

struct haplotypeCohort{
private:
  const linearReferenceStructure* reference;
  size_t number_of_haplotypes;
  bool finalized = false;
  // dim 1: haplotype, dim 2: allele value at site
  vector<vector<alleleValue> > haplotype_alleles_by_site_index;
  // dim 1: site, dim 2: allele value
  vector<vector<size_t> > allele_counts_by_site_index;
  vector<vector<vector<size_t> > > haplotype_indices_by_site_and_allele;
public:
  void populate_allele_counts();
  haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
            const linearReferenceStructure* reference);
  haplotypeCohort(const vector<string>& haplotypes, 
            const linearReferenceStructure* reference);
  haplotypeCohort(size_t cohort_size, const linearReferenceStructure* reference);
  ~haplotypeCohort();
  
  // Removes sites where the minor allele frequency is below [double frequency].
  // For multiallelic sites, the minor allele frequency is taken to be the sum
  // of the minor allele frequencies with respect to the most common allele.
  // This constructs both a new (dynamically allocated) haplotypeCohort as well
  // as a new linearReferenceStructure, which is implicitly returned as

  haplotypeCohort* remove_sites_below_frequency(double frequency) const;
  const linearReferenceStructure* get_reference() const;  
  
  void assign_alleles_at_site(size_t i, vector<alleleValue> alleles_at_site);
  
  size_t size() const;
  size_t get_haplotype_count() const;
  alleleValue get_dominant_allele(size_t site) const;
  size_t get_MAC(size_t site) const;
  size_t sum_MACs() const;
  
  const vector<size_t>& get_matches(size_t site_index, alleleValue a) const;
  vector<size_t> get_non_matches(size_t site_index, alleleValue a) const;
  alleleValue allele_at(size_t site_index, size_t haplotype_index) const;
  size_t number_matching(size_t site_index, alleleValue a) const;
  size_t number_not_matching(size_t site_index, alleleValue a) const;
  bool match_is_rare(size_t site_index, alleleValue a) const;
  
  const vector<alleleValue>& get_alleles_at_site(size_t site_index) const;
  
  vector<size_t> get_active_rows(size_t site, alleleValue a) const;
  rowSet get_active_rowSet(size_t site, alleleValue a) const;

  // 1 : successful
  // 0 : locked
  // -1 : out of order
  int add_record(size_t site);
  // 1 : successful
  // 0 : cohort locked
  // -1 : already assigned
  int set_sample_allele(size_t site, size_t sample, alleleValue a);
  
  void set_column(const vector<alleleValue>& values);
  void set_column(const vector<alleleValue>& values, size_t site);
  
  void simulate_read_query(const char* ref_seq,
                           double mutation_rate,
                           double recombination_rate,
                           double uncertainty_rate,
                           size_t* return_read_sites,
                           char* return_read_seq) const;
                           
  void simulate_read_query_2(
                           const char* ref_seq,
                           double mutation_rate,
                           double recombination_rate,
                           double uncertainty_rate,
                           size_t** return_read_sites,
                           size_t* n_read_sites,
                           char* return_read_seq,
                           char** r_s_alleles_1,
                           char** r_s_alleles_2) const; 
};

#endif