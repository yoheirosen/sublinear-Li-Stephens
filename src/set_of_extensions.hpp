#ifndef EXTENSION_SET_H
#define EXTENSION_SET_H

#include "reference.hpp"
#include "DP_map.hpp"
#include "probability.hpp"
#include "penalty_set.hpp"
#include <vector>

using namespace std;

// An extensionSet holds the results of expensive but not haplotype-specific
// calls. Faster lookup that using a rowSet and also holds the current_maps for
// each allele

struct extensionSet{
private:
  vector<DPUpdateMap> current_map;
  vector<rowSet const*> active_rows;
  vector<bool> match_is_rare;
public:
  extensionSet(haplotypeCohort* cohort, size_t site_index);
  
  bool                  get_match_is_rare(size_t i) const;
  alleleValue           get_allele(size_t i) const;
  const rowSet&         get_active_rows(size_t i) const;
  
  void extend_probability_by_allele(fastFwdAlgState* hap_mat, size_t i);
};

#endif