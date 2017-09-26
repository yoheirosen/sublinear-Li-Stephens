#ifndef ROW_SET_H
#define ROW_SET_H

#include <vector>
#include "allele.hpp"

using namespace std;

struct rowSet{
private:
  vector<alleleValue> included_alleles;
  vector<size_t> boundaries = {0};
  vector<const vector<size_t>* > row_vectors;
public:
  rowSet(vector<const vector<size_t>* > row_vectors, vector<alleleValue> allele);
  size_t operator[](size_t i) const;
  size_t size() const;
};

#endif