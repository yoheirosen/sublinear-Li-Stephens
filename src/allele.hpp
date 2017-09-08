#include <string>

using namespace std;

typedef enum alleleValue{
  A,
  C,
  T,
  G,
  gap
} alleleValue;

char allele_to_char(alleleValue a);
// converts unexpected input to ref
alleleValue char_to_allele(char c, alleleValue ref);
// does not handle unexpected input
alleleValue char_to_allele(char c);
alleleValue str_to_allele(string& s);

struct alleleAtSite{
  size_t site_index;
  alleleValue allele;
  alleleAtSite(size_t site, alleleValue allele);
  alleleAtSite();
};