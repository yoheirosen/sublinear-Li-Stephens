#ifndef LINEAR_INPUT_HAPLOTYPE_H
#define LINEAR_INPUT_HAPLOTYPE_H

#include <string>
#include <vector>
#include <iostream>
#include "reference.hpp"

using std::vector;

struct inputHaplotype{
private:
  siteIndex *reference = NULL;
  // Alleles for each site along the linear reference that this path touches
  // Site identity is inferred from haplotype start and end positions
  vector<alleleValue> alleles;
  // For each inter-allele gap, stores the number of novel SNVs that occur in the haplotype relative to the reference.
  // Includes the leading and trailing gaps.
  vector<size_t> novel_SNVs;
  
  // absolute positions of start and end
  size_t absolute_start_pos;
  size_t absolute_end_pos;
  size_t start_offset_wrt_ref;
  // indices w.r.t. the reference structure of the first and last sites covered
  // by the haplotype (the haplotype may begin and end in the middle of a span)
  size_t start_site = 0;
  size_t end_site = 0;
  // first and last span lengths
  size_t left_tail_length;
  size_t right_tail_length;
  // since we index things by site; we need a flag for the non-existence of any
  // sites within the haplotype.
  // If false, alleles must not be empty.
  bool has_no_sites = false;
  bool invalid = false;
  
  void build(const char* query, const char* reference_sequence, size_t length);
  void calculate_relative_positions(bool covers_reference);
  
  // binary search for site coming before position p. In order to give a
  // meaningful answer, there must be a site below p
  size_t find_site_below(size_t p) const;
public:
  inputHaplotype();
  inputHaplotype(siteIndex* reference);
  inputHaplotype(const vector<alleleValue>& query);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count);
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count,
            siteIndex *reference);          
  inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count,
            siteIndex *reference, size_t absolute_start_pos, 
            size_t length);
  inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference);          
  inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference, size_t absolute_start_pos, 
            size_t length);
  ~inputHaplotype();

#ifdef DEBUG
  void print(ostream& out) const;
  size_t get_length() const;
  bool is_valid() const;
  void validate() const;
#endif            
  
  alleleValue get_allele(size_t j) const;
  const vector<alleleValue>& get_alleles() const;
  size_t get_start_site() const;
  
  // Get the number of novel SNVs in the given gap between alleles (where -1 is the space before the first allele)
  size_t get_n_novel_SNVs(int j) const;
  
  // Length of the haplotype before encountering the first allele
  size_t get_left_tail() const;
  bool has_left_tail() const;
  // Get the distance from this allele to the next, or to the end of the read if it is the last allele
  size_t get_span_after(size_t i) const;
  // True if there is a next allele or span from the allele's end to the end of the read after the allele at the given site.
  bool has_span_after(size_t i) const;
  
  size_t get_site_index(size_t j) const;
  bool has_sites() const;
  size_t number_of_sites() const;
};

#endif