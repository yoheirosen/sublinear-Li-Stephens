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
  unordered_map<size_t, size_t> position_to_site_index;
  vector<size_t> site_index_to_position;
  vector<alleleValue> site_index_to_reference_allele;
  
  void add_site(size_t position, alleleValue reference_value);
  void calculate_final_span_length(size_t reference_length);
public:
  linearReferenceStructure(const vector<string>& haplotypes,
            const string& reference_values);
  linearReferenceStructure(const vector<size_t>& positions, size_t length,
            const vector<alleleValue>& reference_values);
  linearReferenceStructure(const vector<size_t>& positions,
            const string& reference_sequence);
  ~linearReferenceStructure();
  
  vector<size_t> span_lengths;
  size_t leading_span_length;
  
  bool is_site(size_t actual_position) const;
  size_t get_site_index(size_t actual_position) const;
  size_t get_position(size_t site_index) const;
  
  bool has_span_before(size_t site_index) const;
  bool has_span_after(size_t site_index) const;
  size_t span_length_before(size_t site_index) const;
  size_t span_length_after(size_t site_index) const;
  
  bool is_augmentation(alleleValue a, size_t position) const;
  
  size_t number_of_sites() const;
  size_t absolute_length() const;
  
  alleleValue get_reference_allele_at_site(size_t site_index) const;
  
  // behavior when there is no site above: returns the 1-past-the-end index
  size_t find_site_above(size_t position) const;
  // behavior when there is no site below: returns SIZE_MAX
  size_t find_site_below(size_t position) const;
};

struct haplotypeCohort{
private:
  const linearReferenceStructure* reference;
  size_t number_of_haplotypes;
  // dim 1: haplotype, dim 2: allele value at site
  vector<vector<alleleValue> > haplotype_alleles_by_site_index;
  // dim 1: site, dim 2: allele value
  vector<vector<size_t> > allele_counts_by_site_index;
  vector<vector<vector<size_t> > > haplotype_indices_by_site_and_allele;
  void populate_allele_counts();
public:
  haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
            const linearReferenceStructure* reference);
  haplotypeCohort(const vector<string>& haplotypes, 
            const linearReferenceStructure* reference);
  haplotypeCohort(size_t cohort_size, const linearReferenceStructure* reference);
  ~haplotypeCohort();
  
  void assign_alleles_at_site(size_t i, vector<alleleValue> alleles_at_site);
  
  size_t size() const;
  
  const vector<size_t>& get_matches(size_t site_index, alleleValue a) const;
  vector<size_t> get_non_matches(size_t site_index, alleleValue a) const;
  alleleValue allele_at(size_t site_index, size_t haplotype_index) const;
  size_t number_matching(size_t site_index, alleleValue a) const;
  size_t number_not_matching(size_t site_index, alleleValue a) const;
  bool match_is_rare(size_t site_index, alleleValue a) const;
  
  const vector<alleleValue>& get_alleles_at_site(size_t site_index) const;
  
  vector<size_t> get_active_rows(size_t site, alleleValue a) const;
  rowSet get_active_rowSet(size_t site, alleleValue a) const;
};

#endif