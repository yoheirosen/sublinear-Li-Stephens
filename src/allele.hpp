#ifndef ALLELE_H
#define ALLELE_H

#include <string>

using namespace std;

const char N_VALID_ALLELES = 5;

typedef enum alleleValue{
  A,
  C,
  T,
  G,
  gap,
  unassigned
} alleleValue;

namespace allele{
  char to_char(alleleValue a);
  // converts unexpected input to ref
  alleleValue from_char(char c, alleleValue ref);
  // does not handle unexpected input
  alleleValue from_char(char c);
}

struct alleleAtSite{
  size_t site_index;
  alleleValue allele;
  alleleAtSite(size_t site, alleleValue allele);
  alleleAtSite();
};

#endif