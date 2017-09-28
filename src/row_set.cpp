#include "row_set.hpp"

rowSet::rowSet() {
  
}

rowSet::rowSet(vector<const vector<size_t>* > row_vectors, vector<alleleValue> alleles) : included_alleles(alleles), row_vectors(row_vectors) {
  for(size_t i = 1; i < row_vectors.size(); i++) {
    boundaries.push_back(boundaries.back() + row_vectors[i - 1]->size());
  }
}

size_t rowSet::operator[](size_t i) const {
  for(size_t j = boundaries.size() - 1; j >= 0; j--) {
    if(i >= boundaries[j]) {
      i -= boundaries[j];
      return row_vectors[j]->at(i);
    }
  }
}

size_t rowSet::size() const {
  return boundaries.back() + row_vectors.back()->size();
}