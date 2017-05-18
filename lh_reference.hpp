#ifndef LINEAR_REFERENCE_STRUCTURE_H
#define LINEAR_REFERENCE_STRUCTURE_H

#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

typedef enum alleleValue{
  A,
  C,
  T,
  G,
  gap
} alleleValue;

struct linearReferenceStructure{
private:
  unordered_map<size_t, size_t> position_to_site_index;
  vector<size_t> site_index_to_position;
  vector<alleleValue> site_index_to_reference_allele;
  
  void add_site(size_t position, alleleValue reference_value);
  void calculate_final_span_length(size_t reference_length);
public:
  linearReferenceStructure(vector<string>& haplotypes,
            string& reference_values);
  linearReferenceStructure(vector<size_t>& positions, size_t length,
            vector<alleleValue>& reference_values);
  ~linearReferenceStructure();
  
  vector<size_t> span_lengths;
  size_t leading_span_length;
  
  bool is_site(size_t actual_position);
  size_t get_site_index(size_t actual_position);
  size_t get_position(size_t site_index);
  
  bool has_span_before(size_t site_index);
  bool has_span_after(size_t site_index);
  size_t span_length_before(size_t site_index);
  size_t span_length_after(size_t site_index);
  
  bool is_augmentation(alleleValue a, size_t position);
  
  size_t number_of_sites();
  size_t absolute_length();
  
  alleleValue get_reference_allele_at_site(size_t site_index);
};

// converts unexpected input to ref
alleleValue char_to_allele(char c, alleleValue ref);
// does not handle unexpected input
alleleValue char_to_allele(char c);

struct haplotypeCohort{
private:
  linearReferenceStructure* reference;
  size_t number_of_haplotypes;
  vector<vector<alleleValue> > haplotype_alleles_by_site_index;
  vector<vector<size_t> > allele_counts_by_site_index;
  void populate_allele_counts();
public:
  haplotypeCohort(vector<vector<alleleValue> >& haplotypes,
            linearReferenceStructure* reference);
  haplotypeCohort(vector<string>& haplotypes, 
            linearReferenceStructure* reference);
  ~haplotypeCohort();
  
  size_t size();
  
  vector<size_t> get_matches(size_t site_index, alleleValue a);
  alleleValue allele_at(size_t site_index, size_t haplotype_index);
  size_t number_matching(size_t site_index, alleleValue a);
  size_t number_not_matching(size_t site_index, alleleValue a);
};

#endif