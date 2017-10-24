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
    case unassigned:
      return 'N';
  }
}

alleleValue char_to_allele(char c, alleleValue ref) {
  switch (c) {
    case 'a':
    case 'A':
      return A;
    case 'c':
    case 'C':
      return C;
    case 't':
    case 'T':
      return T;
    case 'g':
    case 'G':
      return G;
    case '-':
      return gap;
    default:
      return ref;
  }
}

alleleValue char_to_allele(char c) {
  return char_to_allele(c, unassigned);
}

alleleValue str_to_allele(string& s) {
  return char_to_allele(s[0]);
}


alleleAtSite::alleleAtSite() {
  
}

alleleAtSite::alleleAtSite(size_t site, alleleValue allele) : 
            site_index(site), allele(allele) {
  
}