#ifndef ROW_SET_H
#define ROW_SET_H

#include <vector>
#include "allele.hpp"

using namespace std;

// A rowSet allows us to pass around discontinuous sets of row-indices without
// having to copy subsets of vectors

struct rowSet{
private:
  size_t elements;
  size_t unique_alleles;
  // the least index, within this struct, of the 1st element of each row_vector
  vector<size_t> boundaries = {0};
  vector<const vector<size_t>* > row_vectors;
  vector<alleleValue> included_alleles;
public:
  rowSet();
  rowSet(vector<const vector<size_t>* > row_vectors, vector<alleleValue> allele);
  const size_t& operator[](size_t i) const;
  const size_t& size() const;
};

#endif