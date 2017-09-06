#ifndef LINEAR_INPUT_HAPLOTYPE_H
#define LINEAR_INPUT_HAPLOTYPE_H

#include <string>
#include <vector>
#include "reference.hpp"

using namespace std;

struct inputHaplotype{
private:
  linearReferenceStructure *reference = NULL;
  vector<alleleValue> alleles;
  vector<size_t> augmentations;
  
  // absolute positions of start and end
  size_t start_position;
  size_t end_position;
  // indices w.r.t. the reference structure of the first and last sites covered
  // by the haplotype (the haplotype may begin and end in the middle of a span)
  size_t start_index = 0;
  size_t end_index = 0;
  // since we index things by site; we need a flag for the non-existence of any
  // sites within the haplotype
  bool has_no_sites = false;
  // first and last span lengths
  size_t left_tail_length;
  size_t right_tail_length;
  
  void build_relative_positions();
public:
  inputHaplotype(linearReferenceStructure* reference);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& augmentation_count);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& augmentation_count,
            linearReferenceStructure *reference, size_t start_pos = 1, 
            size_t length = 0);
            
  inputHaplotype(string query, string reference_sequence, 
            linearReferenceStructure* reference, size_t start_pos = 1, 
            size_t length = 0);
  ~inputHaplotype();
              
  void edit(size_t index, alleleValue a);
  void edit(size_t position, char new_c, char old_c, char ref);
  void edit(size_t start_pos, size_t end_pos, string new_string, 
            string old_string, string ref);
  
  alleleValue get_allele(size_t j) const;
  
  size_t get_augmentations(int j) const;
  
  size_t get_left_tail() const;
  bool has_left_tail() const;
  size_t get_span_after(size_t i) const;
  bool has_span_after(size_t i) const;
  
  // binary search for site coming before position p. In order to give a
  // meaningful answer, there must be a site below p
  size_t find_site_below(size_t p) const;
  
  size_t get_site_index(size_t j) const;
  bool has_sites() const;
  size_t number_of_sites() const;
};

#endif