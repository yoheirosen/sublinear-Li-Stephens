#ifndef ALLELE_H
#define ALLELE_H

#include <vector>

using std::vector;
using std::size_t;

const char N_VALID_ALLELES = 5;

struct siteIndex;

typedef enum alleleValue{
  A,
  C,
  T,
  G,
  gap,
  unassigned
} alleleValue;

struct alleleAtSite{
  size_t site_index;
  alleleValue allele;
  alleleAtSite(size_t site, alleleValue allele);
  alleleAtSite();
};

struct alleleVector{
  vector<alleleValue> entries;
  const siteIndex* base_index = nullptr;
  alleleVector();
  alleleVector(const vector<alleleValue>& entries);
  alleleVector(const vector<alleleValue>& entries, const siteIndex* index);
  size_t size() const;
};

namespace allele{
  char to_char(alleleValue a);
  // converts unexpected input to ref
  alleleValue from_char(char c, alleleValue ref);
  // does not handle unexpected input
  alleleValue from_char(char c);
  alleleVector rebase_down(const alleleVector& input, const siteIndex& new_index);
}

#endif