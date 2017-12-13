#include "allele.hpp"

using namespace std;

char allele::to_char(alleleValue a) {
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

alleleValue allele::from_char(char c, alleleValue ref) {
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

alleleValue allele::from_char(char c) {
  return allele::from_char(c, unassigned);
}

alleleAtSite::alleleAtSite() {
  
}

alleleAtSite::alleleAtSite(size_t site, alleleValue allele) : 
            site_index(site), allele(allele) {
  
}

alleleVector::alleleVector() {};

alleleVector::alleleVector(const vector<alleleValue>& entries) : entries(entries) {
}

alleleVector::alleleVector(const vector<alleleValue>& entries, const siteIndex* index) : entries(entries), base_index(index) {
}

size_t alleleVector::size() const {
  return entries.size();
}

