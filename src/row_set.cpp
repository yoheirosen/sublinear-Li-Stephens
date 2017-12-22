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

size_t rowSet::max_element(const vector<double>& values) const {
  double max = values[row_vectors[0]->at(0)];
  size_t max_index = 0;
  for(size_t i = 0; i < row_vectors.size(); i++) {
    for(size_t j = 0; j < row_vectors[i]->size(); j++) {
      if(values[row_vectors[i]->at(j)] > max) {
        max_index = row_vectors[i]->at(j);
        max = values[max_index];
      }
    }
  }
  return max_index;
}

size_t rowSet::min_element(const vector<double>& values) const {
  double min = values[row_vectors[0]->at(0)];
  size_t min_index = 0;
  for(size_t i = 0; i < row_vectors.size(); i++) {
    for(size_t j = 0; j < row_vectors[i]->size(); j++) {
      if(values[row_vectors[i]->at(j)] < min) {
        min_index = row_vectors[i]->at(j);
        min = values[min_index];
      }
    }
  }
  return min_index;
}

rowSet::const_iterator rowSet::begin() const {
  return const_iterator(row_vectors.begin(), row_vectors[0]->begin(), this);
}

rowSet::const_iterator rowSet::end() const {
  return const_iterator(row_vectors.end(), row_vectors.back()->end(), this);
}

const size_t& rowSet::const_iterator::operator*() const {
  return *(inner_itr);
}

bool rowSet::const_iterator::operator==(const rowSet::const_iterator& other) const {
  return inner_itr == other.inner_itr && outer_itr == other.outer_itr;
}

bool rowSet::const_iterator::operator!=(const rowSet::const_iterator& other) const {
  return !(*this == other);
}

rowSet::const_iterator::const_iterator(outer_itr_t outer_itr, inner_itr_t inner_itr, const rowSet* container) : 
  inner_itr(inner_itr), outer_itr(outer_itr), container(container) {
}

rowSet::const_iterator::const_iterator(const const_iterator& other) : inner_itr(other.inner_itr),
  outer_itr(other.outer_itr), container(other.container) {
}

rowSet::const_iterator& rowSet::const_iterator::operator++() {
  ++inner_itr;
  if(inner_itr == (*outer_itr)->end()) {
    ++outer_itr;
  }
  if(outer_itr != container->row_vectors.end()) {
    inner_itr = (*outer_itr)->begin();
  }
  return *this;
}

rowSet::const_iterator rowSet::const_iterator::operator++(int foo) {
  rowSet::const_iterator to_return(*this);
  ++to_return;
  return to_return;
}