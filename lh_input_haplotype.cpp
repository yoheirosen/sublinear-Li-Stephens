#include "lh_input_haplotype.hpp"

using namespace std;

inputHaplotype::~inputHaplotype() {
  
}

inputHaplotype::inputHaplotype(vector<alleleValue> query) : alleles(query) {
  augmentations = vector<size_t>(query.size(), 0);
}

inputHaplotype::inputHaplotype(vector<alleleValue> query, 
            linearReferenceStructure* reference) : reference(reference),
            alleles(query) {
  augmentations = vector<size_t>(query.size(), 0);
}

inputHaplotype::inputHaplotype(vector<alleleValue> query, 
            vector<size_t> augmentation_count) : alleles(query), 
            augmentations(augmentation_count) {
  
}

inputHaplotype::inputHaplotype(vector<alleleValue> query, 
            vector<size_t> augmentation_count, 
            linearReferenceStructure* reference) : reference(reference),
            alleles(query), augmentations(augmentation_count) {

}

inputHaplotype::inputHaplotype(string query, string reference_sequence, 
            linearReferenceStructure* reference) : reference(reference) {
  if(reference->leading_span_length > 0) {
    size_t counter = 0;
    for(size_t i = 0; i < reference->leading_span_length; i++) {
      if(query[i] != reference_sequence[i]) {
        counter++;
      }
    }
    augmentations.push_back(counter);
  } else {
    augmentations.push_back(0);
  }
  for(size_t i = 0; i < reference->number_of_sites(); i++) {
    size_t pos = reference->get_position(i);
    alleles.push_back(char_to_allele(query[pos], 
                reference->get_reference_allele_at_site(i)));
    if(reference->has_span_after(i)) {
      size_t counter = 0;
      for(size_t j = pos + 1; j <= pos + reference->span_length_after(i); j++) {
        if(query[j] != reference_sequence[j]) {
          counter++;
        }
      }
      augmentations.push_back(counter);
    } else {
      augmentations.push_back(0);
    }
  }
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

size_t inputHaplotype::find_site_below(size_t p) {
  if(reference->is_site(p)) {
    return reference->get_site_index(p);
  } else {
    size_t lower_bound = 0;
    size_t upper_bound = reference->number_of_sites() - 1;
    if(p > reference->get_position(upper_bound)) {
      return upper_bound;
    }
    size_t mid = upper_bound/2;
    while(true) {
      if(reference->get_position(mid) > p) {
        if(mid - lower_bound == 1) {
          return lower_bound;
        }
        upper_bound = mid;
        mid = (upper_bound - lower_bound)/2 + lower_bound;
      } else {
        if(upper_bound - mid == 1) {
          return mid;
        }
        lower_bound = mid;
        mid = (upper_bound - lower_bound)/2 + lower_bound;
      }
    }
  }
}

void inputHaplotype::edit(size_t start_pos, size_t end_pos, string new_string, 
          string old_string, string ref) {
  size_t length = end_pos - start_pos + 1;
  for(size_t i = 0; i < length; i++) {
    edit(start_pos + i, new_string[i], old_string[i], ref[i]);
  }
}

alleleValue inputHaplotype::get_allele(size_t j) {
  return alleles[j];
}

size_t inputHaplotype::get_augmentations(int j) {
  return augmentations[j + 1];
}

// void inputHaplotype::build_relative_positions() {
//   start_index = find_start_index();
//   end_index = find_end_index();
//   left_tail_length = reference->get_position(start_index) - start_position;
//   right_tail_length = end_position - reference->get_position(end_index);
// }

// size_t inputHaplotype::get_global_index(size_t i) {
//   return i - start_index;
// }