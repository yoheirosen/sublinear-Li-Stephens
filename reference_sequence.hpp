#ifndef LINEAR_REFERENCE_SEQUENCE_H
#define LINEAR_REFERENCE_SEQUENCE_H

#include <vector>
#include <string>
#include "lh_reference.hpp"

using namespace std;

struct referenceSequence{
private:
  vector<alleleValue> alleles;
public:
  size_t size() const;
  referenceSequence(const vector<alleleValue>& input);
  referenceSequence(const string& input);
  referenceSequence(const char* input);
  bool matches(size_t i, alleleValue a) const;
  bool matches(size_t i, char a) const;
  alleleValue at(size_t i) const;
};

#endif