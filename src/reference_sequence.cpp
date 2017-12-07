#include "reference_sequence.hpp"

size_t referenceSequence::size() const {
  return alleles.size();
}

referenceSequence::referenceSequence(const vector<alleleValue>& input, size_t offset) : alleles(input), offset(offset) {
  
}

referenceSequence::referenceSequence(const char* input, size_t offset) : offset(offset) {
  for(size_t i = 0; input[i] != '\0'; i++) {
    alleles.push_back(allele::from_char(input[i]));
  }
}

referenceSequence::referenceSequence(const string& input, size_t offset) : offset(offset) {
  for(size_t i = 0; i < input.length(); i++) {
    alleles.push_back(allele::from_char(input[i]));
  }
}

bool referenceSequence::matches(size_t i, alleleValue a) const {
  return alleles[i - offset] == a;
}

bool referenceSequence::matches(size_t i, char a) const {
  return alleles[i - offset] == allele::from_char(a, alleles[i - offset]);
}

alleleValue referenceSequence::at(size_t i) const {
  return alleles[i - offset];
}

vector<size_t> referenceSequence::mismatches(const string& other) const {
  return mismatches(other.c_str(), other.size());
}

vector<size_t> referenceSequence::mismatches(const char* other, size_t str_length) const {
  vector<size_t> to_return;
  for(size_t i = 0; i < str_length; i++) {
    if(!matches(i, other[i])) {
      to_return.push_back(i);
    }
  }
  return to_return;
}