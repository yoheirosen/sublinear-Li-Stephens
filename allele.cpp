#include "allele.hpp"

char allele_to_char(alleleValue a) {
  switch(a) {
    case A:
      return 'A';
    case C:
      return 'C';
    case T:
      return 'T';
    case G:
      return 'G';
    case gap:
      return '-';
  }
}

alleleValue char_to_allele(char c, alleleValue ref) {
  switch (c) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'T':
      return T;
    case 'G':
      return G;
    case '-':
      return gap;
    default:
      return ref;
  }
}

alleleValue char_to_allele(char c) {
  return char_to_allele(c, gap);
}

alleleValue str_to_allele(const string& s) {
  return char_to_allele(s[0]);
}


alleleAtSite::alleleAtSite() {
  
}

alleleAtSite::alleleAtSite(size_t site, alleleValue allele) : 
            site_index(site), allele(allele) {
  
}