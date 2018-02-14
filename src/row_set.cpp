#include "row_set.hpp"
#include <iostream>

rowSet::rowSet() {
  row_vectors = {nullptr};
}

rowSet::rowSet(vector<const vector<size_t>* > row_vectors) : row_vectors(row_vectors), n_row_vectors(row_vectors.size()), sizes(vector<size_t>(row_vectors.size())) {
  for(size_t i = 0; i < n_row_vectors; i++) {
    sizes[i] = row_vectors[i]->size();
  }
  if(row_vectors[0]->size() == 0) {
    n_row_vectors = 0;
  }
}

rowSet::const_iterator rowSet::begin() const {
  return const_iterator(
    row_vectors.data(),
    row_vectors[0]->data(),
    row_vectors.size(),
    row_vectors[0]->size(),
    row_vectors);
}

rowSet::const_iterator rowSet::end() const {
  return const_iterator(
    row_vectors.data() + row_vectors.size(),
    row_vectors.back()->data() + row_vectors.back()->size(),
    0,
    0,
    row_vectors);
}

const size_t& rowSet::const_iterator::operator*() const {
  return *inner_itr;
}

bool rowSet::const_iterator::operator!=(const rowSet::const_iterator& other) const {
  return inner_itr != other.inner_itr || outer_itr != other.outer_itr;
}

bool rowSet::const_iterator::operator==(const rowSet::const_iterator& other) const {
  return !(*this != other);
}

rowSet::const_iterator::const_iterator(const const_iterator& other): 
  outer_itr(other.outer_itr), 
  inner_itr(other.inner_itr),
  inner_counter(other.inner_counter),
  outer_counter(other.outer_counter),
  container(other.container)
  {
}

rowSet::const_iterator::const_iterator(outer_itr_t outer_itr, inner_itr_t inner_itr,   size_t outer_idx, size_t inner_idx, const vector<const vector<size_t>*>& container) : 
  outer_itr(outer_itr), 
  inner_itr(inner_itr),
  inner_counter(inner_idx),
  outer_counter(outer_idx),
  container(&container)
  {
}

rowSet::const_iterator& rowSet::const_iterator::operator++() {
  ++inner_itr;
  --inner_counter;
  if(inner_counter == 0) {
    ++outer_itr;
    --outer_counter;
    if(outer_counter != 0) {
      inner_itr = (*outer_itr)->data();
      inner_counter = (*outer_itr)->size();
    }
  }
  return *this;
}

rowSet::const_iterator rowSet::const_iterator::operator++(int foo) {
  rowSet::const_iterator temp(*this);
  ++(*this);
  return temp;
}

bool rowSet::empty() const {
  return n_row_vectors == 0;
}