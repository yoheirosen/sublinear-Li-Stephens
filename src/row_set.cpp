#include "row_set.hpp"

rowSet::rowSet() {
  
}

rowSet::rowSet(vector<const vector<size_t>* > row_vectors, 
               vector<alleleValue> alleles) :
                        included_alleles(alleles), 
                        row_vectors(row_vectors) {
  unique_alleles = row_vectors.size();
  if(row_vectors.size() == 0) {
    elements = 0;
  } else {
    for(size_t i = 1; i < unique_alleles; i++) {
      boundaries.push_back(boundaries.back() + row_vectors[i - 1]->size());
    }
    elements = boundaries.back() + row_vectors.back()->size();
  }
}

const size_t& rowSet::operator[](size_t i) const {
  // TODO: make optimized accessor for loop
  for(size_t j = unique_alleles - 1; j >= 0; j--) {
    if(i >= boundaries[j]) {
      i -= boundaries[j];
      return row_vectors[j]->at(i);
    }
  }
}

const size_t& rowSet::size() const {
  return elements;
}