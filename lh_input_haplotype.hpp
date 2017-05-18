#ifndef LINEAR_INPUT_HAPLOTYPE_H
#define LINEAR_INPUT_HAPLOTYPE_H

#include <string>
#include <vector>
#include "lh_reference.hpp"

using namespace std;

struct inputHaplotype{
private:
  linearReferenceStructure *reference = NULL;
  vector<alleleValue> alleles;
  vector<size_t> augmentations;
public:
  inputHaplotype(vector<alleleValue> query, vector<size_t> augmentation_count);
  inputHaplotype(vector<alleleValue> query, vector<size_t> augmentation_count,
            linearReferenceStructure *reference);
  inputHaplotype(vector<alleleValue> query);
  inputHaplotype(vector<alleleValue> query,
            linearReferenceStructure *reference);
  inputHaplotype(string query, string reference_sequence, 
              linearReferenceStructure* reference);
  ~inputHaplotype();
              
  void edit(size_t index, alleleValue a);
  void edit(size_t position, char new_c, char old_c, char ref);
  void edit(size_t start_pos, size_t end_pos, string new_string, 
            string old_string, string ref);
  
  bool validate_haplotype();
  bool validate_augmentations();
  alleleValue get_allele(size_t j);
  size_t get_augmentations(int j);
};

#endif