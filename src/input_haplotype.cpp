#include "input_haplotype.hpp"

using namespace std;

inputHaplotype::~inputHaplotype() {
  
}

inputHaplotype::inputHaplotype(linearReferenceStructure* reference) : 
          reference(reference) {
  
}

inputHaplotype::inputHaplotype(const vector<alleleValue>& query, 
            const vector<size_t>& augmentation_count) : alleles(query), 
            augmentations(augmentation_count) {
  
}

inputHaplotype::inputHaplotype(const vector<alleleValue>& query, 
            const vector<size_t>& augmentation_count, 
            linearReferenceStructure* reference, size_t ih_start, 
            size_t length) : reference(reference), alleles(query), 
            augmentations(augmentation_count), ih_start_position(ih_start)
            {
  ih_end_position = ih_start + length - 1;
  build_relative_positions();
}

void inputHaplotype::build_relative_positions() {
  size_t ref_end_site = reference->number_of_sites() - 1; 
  if(ih_end_position == 0 && ih_start_position == 1) {
    // This is impossible and therefore is used as a flag for taking the
    // haplotype to cover the entire length of the reference
    ih_start_position = reference->start_position();
    ih_end_position = reference->start_position() + reference->absolute_length() - 1;
    start_index = 0;
    end_index = ref_end_site;
    left_tail_length = reference->span_length_before(0);
    right_tail_length = reference->span_length_after(ref_end_site);
    return;    
  } else if(!(reference->is_site(ih_start_position)) &&
            !(reference->is_site(ih_end_position))) {
    // Neither the start nor end position are sites. If this is true then it is
    // possible that there are no sites within the interval specified by the
    // inputHaplotype. This can arise in three ways, which we should check
    if(ih_start_position > reference->get_position(ref_end_site) ||
                ih_end_position < reference->get_position(0) ||
                find_site_below(ih_start_position) ==
                find_site_below(ih_end_position)) {
      has_no_sites = true;
      left_tail_length = ih_end_position - ih_start_position + 1;
      return;
    }
  }
  if(reference->is_site(ih_start_position)) {
    start_index = reference->get_site_index(ih_start_position);
    left_tail_length = 0;
  } else {
    // Since we need all indices to have positions within the interval spanned
    // by the haplotype, then a ih_start_position which is not a site must be
    // placed before the 0-index with respect to the haplotype interval 
    if(ih_start_position < reference->get_position(0)) {
      start_index = 0;
      // site positions are not included in span lengths
      left_tail_length = reference->get_position(0) - ih_start_position;
    } else {
     // site_below + 1 is guaranteed to be within range
     start_index = find_site_below(ih_start_position) + 1;
     left_tail_length = reference->get_position(start_index) -
                ih_start_position;
   }
  }
  if(reference->is_site(ih_end_position)) {
    end_index = reference->get_site_index(ih_end_position);
    right_tail_length = 0;
  } else {
    end_index = find_site_below(ih_end_position);
    right_tail_length = ih_end_position - reference->get_position(end_index);
  }
  return;
}

inputHaplotype::inputHaplotype(const char* query, const char* reference_sequence, 
            linearReferenceStructure* reference, size_t ih_start, 
            size_t length) : reference(reference), ih_start_position(ih_start) {
  ih_end_position = ih_start + length - 1;
  ih_start_offset = ih_start_position - reference->start_position();
  build_relative_positions();
  
  size_t number_of_sites;
  if(has_no_sites) {
    number_of_sites = 0;
  } else {
    number_of_sites = end_index - start_index + 1;
  }
  
  if(left_tail_length > 0) {
    size_t counter = 0;
    for(size_t i = 0; i < left_tail_length; i++) {
      if(char_to_allele(query[i]) != char_to_allele(reference_sequence[i + ih_start_offset])) {
        counter++;
      }
    }
    augmentations.push_back(counter);
  } else {
    augmentations.push_back(0);
  }
  if(!has_no_sites) {
    for(size_t i = 0; i < number_of_sites; i++) {
      size_t p_q = 
                reference->get_position(get_site_index(i)) - ih_start_position;
      alleles.push_back(char_to_allele(query[p_q]));
      if(reference->has_span_after(get_site_index(i))) {
        size_t span_end_q;
        if(i == number_of_sites - 1) {
          span_end_q = ih_end_position - ih_start_position;
        } else {
          span_end_q = reference->get_position(get_site_index(i + 1))
                    - 1 - ih_start_position;
        }
        size_t counter = 0;
        for(size_t j = p_q + 1; j <= span_end_q; j++) {
          if(query[j] != reference_sequence[j + ih_start_offset]) {
            counter++;
          }
        }
        augmentations.push_back(counter);
      } else {
        augmentations.push_back(0);
      }
    }
  }
}

bool inputHaplotype::has_sites() const {
  return !has_no_sites;
}

void inputHaplotype::edit(size_t index, alleleValue a) {
  alleles[index] = a;
}

void inputHaplotype::edit(size_t position, char new_c, char old_c, char ref) {
  if(reference->is_site(position)) {
    alleles[reference->get_site_index(position)] = 
              char_to_allele(new_c, char_to_allele(ref));
  } else {
    int delta_aug;
    // does the number of augmentations in the span change?
    if(new_c != old_c) {
      if(old_c == ref) {
        delta_aug = 1;
      } else {
        if(new_c == ref) {
          delta_aug = -1;
        } else {
          delta_aug = 0;
        }
      }
    } else {
      delta_aug = 0;
    }
    if(position < reference->get_position(0)) {
      augmentations[0] += delta_aug;
    } else {
      size_t site_index = find_site_below(position);
      augmentations[site_index + 1] += delta_aug;
    }
  }
}

size_t inputHaplotype::find_site_below(size_t p) const {
  return reference->find_site_below(p);
}

void inputHaplotype::edit(size_t start_pos, size_t end_pos, string new_string, 
          string old_string, string ref) {
  size_t length = end_pos - start_pos + 1;
  for(size_t i = 0; i < length; i++) {
    edit(start_pos + i, new_string[i], old_string[i], ref[i]);
  }
}

alleleValue inputHaplotype::get_allele(size_t j) const {
  return alleles[j];
}

size_t inputHaplotype::get_augmentations(int j) const {
  return augmentations[j + 1];
}

size_t inputHaplotype::get_site_index(size_t j) const {
  return j + start_index;
}

size_t inputHaplotype::get_left_tail() const {
  return left_tail_length;
}

bool inputHaplotype::has_left_tail() const {
  return left_tail_length != 0;
}

size_t inputHaplotype::get_span_after(size_t i) const {
  if(i = end_index) {
    return right_tail_length;
  } else {
    return reference->span_length_after(get_site_index(i));
  }
}

bool inputHaplotype::has_span_after(size_t i) const {
  if(i = end_index) {
    return right_tail_length != 0;
  } else {
    return reference->span_length_after(get_site_index(i));
  }
}

size_t inputHaplotype::number_of_sites() const {
  if(has_no_sites) {
    return 0;
  } else {
    return end_index - start_index + 1;
  }
}