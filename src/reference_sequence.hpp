#ifndef LINEAR_REFERENCE_SEQUENCE_H
#define LINEAR_REFERENCE_SEQUENCE_H

#include <vector>
#include <string>
#include "reference.hpp"

using std::vector;

struct referenceSequence{
private:
  vector<alleleValue> alleles;
  size_t offset;
public:
  size_t size() const;
  referenceSequence(const vector<alleleValue>& input, size_t offset = 0);
  referenceSequence(const string& input, size_t offset = 0);
  referenceSequence(const char* input, size_t offset = 0);
  bool matches(size_t i, alleleValue a) const;
  bool matches(size_t i, char a) const;
  alleleValue at(size_t i) const;
  vector<size_t> mismatches(const string& other) const;
  vector<size_t> mismatches(const char* other, size_t str_length) const;
};

#endif