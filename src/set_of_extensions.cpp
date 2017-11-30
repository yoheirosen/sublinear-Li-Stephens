#include "set_of_extensions.hpp"

extensionSet::extensionSet(haplotypeCohort* cohort, size_t site_index) {
  alleleValue a;
  for(size_t i = 0; i < 5; i++) {
    a = get_allele(i);
    match_is_rare.push_back(cohort->match_is_rare(site_index, a));
    active_rows.push_back(cohort->get_active_rows(site_index, a));
  }
}

bool extensionSet::get_match_is_rare(size_t i) const {
  return match_is_rare[i];
}

alleleValue extensionSet::get_allele(size_t i) const  {
  return (alleleValue)i;
}

const vector<size_t>& extensionSet::get_active_rows(size_t i) const  {
  return active_rows[i];
}

void extensionSet::extend_probability_by_allele(fastFwdAlgState* hap_mat,
            size_t i) {
  hap_mat->extend_probability_at_site(current_map[i], active_rows[i],
              match_is_rare[i], get_allele(i));
}
