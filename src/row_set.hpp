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
  rowSet(const vector<size_t>* row_vector, alleleValue allele);
  const size_t& operator[](size_t i) const;
  const size_t& size() const;

  struct const_iterator{
  public:
    typedef vector<size_t>::const_iterator inner_itr_t;
    typedef vector<const vector<size_t>* >::const_iterator outer_itr_t;
    const_iterator(outer_itr_t outer_itr, inner_itr_t inner_itr, const rowSet* container);
    const_iterator(const const_iterator& other);
    const_iterator& operator++();
    const_iterator operator++(int foo);
    const size_t& operator*() const;
    bool operator==(const rowSet::const_iterator& other) const;
    bool operator!=(const rowSet::const_iterator& other) const;
  private:
    inner_itr_t inner_itr;
    outer_itr_t outer_itr; 
    const rowSet* container;
  };

  const_iterator begin() const;
  const_iterator end() const;
};

#endif