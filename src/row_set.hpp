
#ifndef ROW_SET_H
#define ROW_SET_H

#include <vector>
#include "allele.hpp"

using namespace std;

// A rowSet allows us to pass around discontinuous sets of row-indices without
// having to copy subsets of vectors

struct rowSet{
private:
  size_t n_row_vectors = 0;
  // the least index, within this struct, of the 1st element of each row_vector
  vector<size_t> sizes = {0};
  vector<const vector<size_t>* > row_vectors;
public:
  rowSet();
  rowSet(vector<const vector<size_t>* > row_vectors);

  struct const_iterator{
  public:
    // typedef size_t inner_itr_t;
    // typedef size_t outer_itr_t;
    typedef const size_t* inner_itr_t;
    typedef const vector<size_t>* const* outer_itr_t;
    const_iterator(outer_itr_t outer_itr, inner_itr_t inner_itr, size_t outer_idx, size_t inner_idx, const vector<const vector<size_t>*>& container);
    const_iterator(const const_iterator& other);
    const_iterator& operator++();
    const_iterator operator++(int foo);
    const size_t& operator*() const;
    bool operator==(const rowSet::const_iterator& other) const;
    bool operator!=(const rowSet::const_iterator& other) const;
  private:
    size_t inner_counter;
    size_t outer_counter;
    inner_itr_t inner_itr;
    outer_itr_t outer_itr;
    // const rowSet* container;
    const vector<const vector<size_t>* >* container;
  };

  const_iterator begin() const;
  const_iterator end() const;
  bool empty() const;
};

#endif