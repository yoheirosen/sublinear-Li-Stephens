#ifndef LINEAR_INPUT_HAPLOTYPE_H
#define LINEAR_INPUT_HAPLOTYPE_H

#include <string>
#include <vector>
#include "reference.hpp"

using namespace std;

struct inputHaplotype{
private:
  siteIndex *reference = NULL;
  vector<alleleValue> alleles;
  vector<size_t> augmentations;
  
  // absolute positions of start and end
  size_t ih_start_position;
  size_t ih_end_position;
  size_t ih_start_offset;
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
  
  void build(const char* query, const char* reference_sequence, size_t length);
  void build_relative_positions();
  
  // binary search for site coming before position p. In order to give a
  // meaningful answer, there must be a site below p
  size_t find_site_below(size_t p) const;
public:
  inputHaplotype();
  inputHaplotype(siteIndex* reference);
  inputHaplotype(const vector<alleleValue>& query);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& augmentation_count);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& augmentation_count,
            siteIndex *reference);          
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& augmentation_count,
            siteIndex *reference, size_t ih_start_position, 
            size_t length);
  inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference);          
  inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference, size_t ih_start_position, 
            size_t length);
  ~inputHaplotype();
              
  alleleValue get_allele(size_t j) const;
  const vector<alleleValue>& get_alleles() const;
  size_t get_start_index() const;
  
  size_t get_augmentations(int j) const;
  
  size_t get_left_tail() const;
  bool has_left_tail() const;
  size_t get_span_after(size_t i) const;
  bool has_span_after(size_t i) const;
  
  size_t get_site_index(size_t j) const;
  bool has_sites() const;
  size_t number_of_sites() const;
};

#endif