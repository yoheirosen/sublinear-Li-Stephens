#include "input_haplotype.hpp"
#include <iostream>
#include <cassert>

using namespace std;

inputHaplotype::~inputHaplotype() {
  
}

inputHaplotype::inputHaplotype(siteIndex* reference) : 
          reference(reference) {
  
}

inputHaplotype::inputHaplotype() : has_no_sites(true), invalid(true) {
  
}

inputHaplotype::inputHaplotype(const vector<alleleValue>& query) : alleles(query), 
            novel_SNVs(vector<size_t>(query.size(), 0)) {
  
}

inputHaplotype::inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count) : 
            alleles(query), novel_SNVs(novel_SNV_count) {
  
}

inputHaplotype::inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count, 
            siteIndex* reference, size_t absolute_start_pos, size_t length = 0) : reference(reference), alleles(query), 
            novel_SNVs(novel_SNV_count), absolute_start_pos(absolute_start_pos)
            {
  absolute_end_pos = absolute_start_pos + length - 1;
  calculate_relative_positions(false);
  // Make sure we got a consistent answer about whether we have sites or not.
  assert(has_no_sites == (alleles.size() == 0));
}

//TODO check if ever called
inputHaplotype::inputHaplotype(const vector<alleleValue>& query, const vector<size_t>& novel_SNV_count, 
            siteIndex* reference) : reference(reference), alleles(query), novel_SNVs(novel_SNV_count) {
  calculate_relative_positions(true);
}

//TODO: bounds checking wrt reference start and finish
void inputHaplotype::calculate_relative_positions(bool covers_reference) {
  size_t first_ref_site = 0;
  size_t last_ref_site = reference->number_of_sites() - 1;
  size_t first_ref_site_pos = reference->get_position(0);
  size_t last_ref_site_pos = reference->get_position(last_ref_site);
  size_t ref_start_pos = reference->start_position();
  size_t ref_end_pos = reference->end_position();
  
  if(covers_reference) {
    absolute_start_pos = ref_start_pos;
    absolute_end_pos = ref_end_pos;
    start_site = first_ref_site;
    end_site = last_ref_site;
    left_tail_length = reference->span_length_before(0);
    right_tail_length = reference->span_length_after(last_ref_site);
    return;    
  } 
  
  // Check whether the input_haplotype contains zero sites
  if((!(reference->is_site(absolute_start_pos)) && !(reference->is_site(absolute_end_pos))) &&
     (absolute_start_pos > reference->get_position(last_ref_site) ||
      absolute_end_pos < reference->get_position(0) ||
      find_site_below(absolute_start_pos) == find_site_below(absolute_end_pos))) {
    has_no_sites = true;
    left_tail_length = absolute_end_pos - absolute_start_pos + 1;
    right_tail_length = 0;
    return;
  }
  
  if(absolute_start_pos < ref_start_pos ||
     absolute_start_pos > ref_end_pos ||
     absolute_end_pos < ref_start_pos ||
     absolute_end_pos > ref_end_pos) {
    throw runtime_error("haplotype queried is outside of reference range");
  }
  if(absolute_end_pos < absolute_start_pos) {
    throw runtime_error("start position of haplotype queried exceeds end position");
  }
  
  if(reference->is_site(absolute_start_pos)) {
    start_site = reference->get_site_index(absolute_start_pos);
    left_tail_length = 0;
  } else {
    // Since we need all indices to have positions within the interval spanned
    // by the haplotype, then a absolute_start_pos which is not a site must be
    // placed before the 0-index with respect to the haplotype interval 
    if(absolute_start_pos < reference->get_position(0)) {
      start_site = 0;
      // site positions are not included in span lengths
      left_tail_length = reference->get_position(0) - absolute_start_pos;
    } else {
      // site_below + 1 is guaranteed to be within range
      start_site = find_site_below(absolute_start_pos) + 1;
      left_tail_length = reference->get_position(start_site) - absolute_start_pos;
   }
  }
  
  if(reference->is_site(absolute_end_pos)) {
    end_site = reference->get_site_index(absolute_end_pos);
    right_tail_length = 0;
  } else {
    end_site = find_site_below(absolute_end_pos);
    right_tail_length = absolute_end_pos - reference->get_position(end_site);
  }
  
  start_offset_wrt_ref = absolute_start_pos - reference->start_position();
  return;
}

size_t inputHaplotype::find_site_below(size_t p) const {
  return reference->find_site_below(p);
}

void inputHaplotype::build(const char* query, const char* reference_sequence, size_t length) {
  absolute_end_pos = absolute_start_pos + length - 1;
  calculate_relative_positions(false);
  
  size_t number_of_sites;
  if(has_no_sites) {
    number_of_sites = 0;
  } else {
    number_of_sites = end_site - start_site + 1;
  }
  
  if(left_tail_length > 0) {
    size_t counter = 0;
    for(size_t i = 0; i < left_tail_length; i++) {
      if(allele::from_char(query[i]) != allele::from_char(reference_sequence[i + start_offset_wrt_ref])) {
        counter++;
      }
    }
    novel_SNVs.push_back(counter);
  } else {
    novel_SNVs.push_back(0);
  }
  if(!has_no_sites) {
    for(size_t i = 0; i < number_of_sites; i++) {
      size_t p_q = reference->get_position(get_site_index(i)) - absolute_start_pos;
      alleles.push_back(allele::from_char(query[p_q]));
      if(reference->has_span_after(get_site_index(i))) {
        size_t span_end_q;
        if(i == number_of_sites - 1) {
          span_end_q = absolute_end_pos - absolute_start_pos;
        } else {
          span_end_q = reference->get_position(get_site_index(i + 1))
                    - 1 - absolute_start_pos;
        }
        size_t counter = 0;
        for(size_t j = p_q + 1; j <= span_end_q; j++) {
          if(query[j] != reference_sequence[j + start_offset_wrt_ref]) {
            counter++;
          }
        }
        novel_SNVs.push_back(counter);
      } else {
        novel_SNVs.push_back(0);
      }
    }
  }
}

inputHaplotype::inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference, size_t ih_start = 1, 
            size_t length = 0) : reference(reference), absolute_start_pos(ih_start) {
  build(query, reference_sequence, length);
}

inputHaplotype::inputHaplotype(const char* query, const char* reference_sequence, 
            siteIndex* reference) : reference(reference), absolute_start_pos(reference->start_position()) {
  build(query, reference_sequence, reference->length_in_bp());
}

bool inputHaplotype::has_sites() const {
  return !has_no_sites;
}

alleleValue inputHaplotype::get_allele(size_t j) const {
  #ifdef DEBUG
    return alleles.at(j);
  #else
    return alleles[j];
  #endif
}

size_t inputHaplotype::get_n_novel_SNVs(int j) const {
  #ifdef DEBUG
    return novel_SNVs.at(j + 1);
  #else
    return novel_SNVs[j + 1];
  #endif
}

size_t inputHaplotype::get_site_index(size_t j) const {
  return j + start_site;
}

size_t inputHaplotype::get_left_tail() const {
  return left_tail_length;
}

bool inputHaplotype::has_left_tail() const {
  return left_tail_length != 0;
}

size_t inputHaplotype::get_span_after(size_t i) const {
  if(i == end_site) {
    return right_tail_length;
  } else {
    return reference->span_length_after(get_site_index(i));
  }
}

bool inputHaplotype::has_span_after(size_t i) const {
  if(i == end_site) {
    return right_tail_length != 0;
  } else {
    return reference->span_length_after(get_site_index(i));
  }
}

size_t inputHaplotype::number_of_sites() const {
  if(has_no_sites) {
    return 0;
  } else {
    return end_site - start_site + 1;
  }
}

const vector<alleleValue>& inputHaplotype::get_alleles() const {
  return alleles;
}

bool inputHaplotype::is_valid() const {
  return !invalid;
}

void inputHaplotype::validate() const {
  for(size_t i = start_site; i < end_site; i++) {
    alleleValue test = alleles[i - start_site];
    if(test > 5) {
      cout << i << " " << test << endl; 
      throw runtime_error("invalid allele in observed haplotype");
    }
  }
}

size_t inputHaplotype::get_start_site() const {
  return start_site;
}

size_t inputHaplotype::get_length() const {
  return absolute_end_pos - absolute_start_pos + 1;
}