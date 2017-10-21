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
  size_t allele_count;
  // the least index, within this struct, of the 1st element of each row_vector
  vector<size_t> lower_bounds = {0};
  vector<const vector<size_t>* > row_vectors;
  vector<alleleValue> included_alleles;
public:
  rowSet();
  rowSet(vector<const vector<size_t>* > row_vectors, vector<alleleValue> allele);
  const size_t& operator[](size_t i) const;
  const size_t& size() const;

  size_t max_element(const vector<double>& values) const;
  size_t min_element(const vector<double>& values) const;
};

#endif