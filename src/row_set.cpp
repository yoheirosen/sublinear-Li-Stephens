#include "row_set.hpp"

rowSet::rowSet() {
  row_vectors = {nullptr};
}

rowSet::rowSet(vector<const vector<size_t>* > row_vectors, 
               vector<alleleValue> alleles) :
                        included_alleles(alleles), 
                        row_vectors(row_vectors) {
  allele_count = row_vectors.size();
  if(row_vectors.size() == 0) {
    elements = 0;
  } else {
    for(size_t i = 1; i < allele_count; i++) {
      lower_bounds.push_back(lower_bounds.back() + row_vectors[i - 1]->size());
    }
    elements = lower_bounds.back() + row_vectors.back()->size();
  }
}

const size_t& rowSet::operator[](size_t i) const {
  // TODO (10/3/17): make optimized accessor for loop
  for(size_t j = allele_count - 1; j >= 0; j--) {
    if(i >= lower_bounds[j]) {
      i -= lower_bounds[j];
      return row_vectors[j]->at(i);
    }
  }
}

const size_t& rowSet::size() const {
  return elements;
}

