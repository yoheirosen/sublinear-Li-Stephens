#include "reference_sequence.hpp"

size_t referenceSequence::size() const {
  return alleles.size();
}

referenceSequence::referenceSequence(const vector<alleleValue>& input) : alleles(input) {
  
}

referenceSequence::referenceSequence(const char* input) {
  for(size_t i = 0; input[i] != '\0'; i++) {
    alleles.push_back(char_to_allele(input[i]));
  }
}

referenceSequence::referenceSequence(const string& input) {
  for(size_t i = 0; i < input.length(); i++) {
    alleles.push_back(char_to_allele(input[i]));
  }
}

bool referenceSequence::matches(size_t i, alleleValue a) const {
  return alleles[i] == a;
}

bool referenceSequence::matches(size_t i, char a) const {
  return alleles[i] == char_to_allele(a, alleles[i]);
}

alleleValue referenceSequence::at(size_t i) const {
  return alleles[i];
}