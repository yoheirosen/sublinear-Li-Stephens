#include "lh_reference.hpp"

using namespace std;

linearReferenceStructure::~linearReferenceStructure() {
  
}

bool linearReferenceStructure::is_site(size_t actual_position) const {
  return (position_to_site_index.count(actual_position) == 1);
}

size_t linearReferenceStructure::get_site_index(size_t actual_position) const {
  return position_to_site_index.at(actual_position);
}

size_t linearReferenceStructure::get_position(size_t site_index) const {
  return site_index_to_position[site_index];
}

bool linearReferenceStructure::has_span_before(size_t site_index) const {
  if(site_index == 0) {
    return (leading_span_length != 0);
  } else {
    return (span_lengths[site_index - 1] != 0);
  }
}

bool linearReferenceStructure::has_span_after(size_t site_index) const {
  return (span_lengths[site_index] != 0);
}

size_t linearReferenceStructure::span_length_before(size_t site_index) const {
  if(site_index == 0) {
    return leading_span_length;
  } else {
    return span_lengths[site_index - 1];
  }
}

size_t linearReferenceStructure::span_length_after(size_t site_index) const {
  return span_lengths[site_index];
}

size_t linearReferenceStructure::number_of_sites() const {
  return site_index_to_position.size();
}

size_t linearReferenceStructure::absolute_length() const {
  return site_index_to_position.back() + span_lengths.back() + 1;
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

linearReferenceStructure::linearReferenceStructure(
            const vector<string>& haplotypes,
            const string& reference_values) {
  size_t number_of_haplotypes = haplotypes.size();
  size_t length_of_haplotypes = haplotypes[0].size();
  for(size_t i = 0; i < length_of_haplotypes; i++) {
    vector<alleleValue> alleles_seen;
    alleleValue ref_at_i = char_to_allele(reference_values[i], A);
    alleles_seen.push_back(ref_at_i);
    for(size_t j = 0; j < number_of_haplotypes; j++) {
      alleleValue allele_at_j = char_to_allele(haplotypes[j][i], ref_at_i);
      bool allele_seen = false;
      for(size_t k = 0; k < alleles_seen.size(); k++) {
        if(allele_at_j == alleles_seen[k]) {
          allele_seen = true;
        }
      }
      if(!allele_seen) {
        alleles_seen.push_back(allele_at_j);
      }
    }
    if(alleles_seen.size() > 1) {
      add_site(i, ref_at_i);
    }
  }
  calculate_final_span_length(length_of_haplotypes);
}

linearReferenceStructure::linearReferenceStructure(
            const vector<size_t>& positions,
            size_t length, const vector<alleleValue>& reference_values) {
  for(size_t i = 0; i < positions.size(); i++) {
    add_site(positions[i], reference_values[i]);
  }
  calculate_final_span_length(length);
}

linearReferenceStructure::linearReferenceStructure(
            const vector<size_t>& positions,
            const string& reference_sequence) {
  alleleValue allele;
  for(size_t i = 0; i < positions.size(); i++) {
    allele = char_to_allele(reference_sequence[positions[i]]);
    add_site(positions[i], allele);
  }
  calculate_final_span_length(reference_sequence.length());
}


void linearReferenceStructure::add_site(size_t position, 
            alleleValue reference_value) {
  size_t new_index = site_index_to_position.size();
  if(new_index == 0) {
    leading_span_length = position;
  } else {
    size_t previous_position = site_index_to_position.back();
    span_lengths.push_back(position - previous_position - 1);
  }
  site_index_to_position.push_back(position);
  site_index_to_reference_allele.push_back(reference_value);
  position_to_site_index.emplace(position, new_index);
}

void linearReferenceStructure::calculate_final_span_length(
            size_t reference_length) {
  size_t previous_position = site_index_to_position.back();
  span_lengths.push_back(reference_length - previous_position - 1);
}

alleleValue linearReferenceStructure::get_reference_allele_at_site(size_t site_index) const {
  return site_index_to_reference_allele[site_index];
}

size_t linearReferenceStructure::find_site_above(size_t position) const {
  if(is_site(position)) {
    return get_site_index(position);
  } else {
    if(get_position(0) > position) {
      return 0;
    } else {
      return find_site_below(position) + 1;
    }
  }
}

size_t linearReferenceStructure::find_site_below(size_t position) const {
  if(is_site(position)) {
    return get_site_index(position);
  } else {
    size_t lower_bound = 0;
    size_t upper_bound = number_of_sites() - 1;
    if(position < get_position(0)) {
      return SIZE_MAX;
    }
    if(position > get_position(upper_bound)) {
      return upper_bound;
    }
    size_t mid = upper_bound/2;
    while(true) {
      if(get_position(mid) > position) {
        if(mid - lower_bound == 1) {
          return lower_bound;
        }
        upper_bound = mid;
        mid = (upper_bound - lower_bound)/2 + lower_bound;
      } else {
        if(upper_bound - mid == 1) {
          return mid;
        }
        lower_bound = mid;
        mid = (upper_bound - lower_bound)/2 + lower_bound;
      }
    }
  }
}

haplotypeCohort::~haplotypeCohort() {
  
}

alleleValue haplotypeCohort::allele_at(size_t site_index, size_t haplotype_index) const {
  return haplotype_alleles_by_site_index[haplotype_index][site_index];
}

size_t haplotypeCohort::number_matching(size_t site_index, alleleValue a) const {
  return allele_counts_by_site_index[site_index][(int)a];
}

size_t haplotypeCohort::number_not_matching(size_t site_index, alleleValue a) const {
  return number_of_haplotypes - number_matching(site_index, a);
}

void haplotypeCohort::populate_allele_counts() {
  int num_sites = reference->number_of_sites();
  allele_counts_by_site_index = vector<vector<size_t> >(num_sites,
              vector<size_t>(5, 0));
  haplotype_indices_by_site_and_allele = vector<vector<vector<size_t> > >(
              num_sites, vector<vector<size_t> >(5, vector<size_t>()));
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    for(int j = 0; j < num_sites; j++) {
      int allele_rank = (int)haplotype_alleles_by_site_index[i][j];
      allele_counts_by_site_index[j][allele_rank]++;
      haplotype_indices_by_site_and_allele[j][allele_rank].push_back(i);
    }
  }
}

size_t haplotypeCohort::size() const {
  return number_of_haplotypes;
}

const vector<size_t>& haplotypeCohort::get_matches(size_t site_index, alleleValue a) const {
  size_t allele_rank = (size_t)a;
  return haplotype_indices_by_site_and_allele[site_index][allele_rank];
}

haplotypeCohort::haplotypeCohort(const vector<string>& haplotypes, 
            const linearReferenceStructure* reference) : reference(reference) {
  number_of_haplotypes = haplotypes.size();
  size_t num_sites = reference->number_of_sites();
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    haplotype_alleles_by_site_index.push_back(vector<alleleValue>(num_sites));
    for(size_t j = 0; j < num_sites; j++) {
      size_t pos = reference->get_position(j);
      haplotype_alleles_by_site_index[i][j] = char_to_allele(haplotypes[i][pos],
                  reference->get_reference_allele_at_site(j));
    }
  }
  populate_allele_counts();
}

haplotypeCohort::haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
          const linearReferenceStructure* reference) : reference(reference),
          haplotype_alleles_by_site_index(haplotypes) {
  number_of_haplotypes = haplotypes.size();
  populate_allele_counts();
}

bool haplotypeCohort::match_is_rare(size_t site_index, alleleValue a) const {
  return number_matching(site_index, a) < number_not_matching(site_index, a);
}

vector<size_t> haplotypeCohort::get_non_matches(size_t site_index, alleleValue a) const {
  vector<size_t> nonmatchlist;
  size_t allele_rank = (size_t)a;
  for(size_t j = 0; j < 5; j++) {
    if(j != allele_rank) {
      for(size_t i = 0; i < 
            (haplotype_indices_by_site_and_allele[site_index][j]).size();
            i++) {
        nonmatchlist.push_back(
              haplotype_indices_by_site_and_allele[site_index][j][i]);
      }
    }
  }
  return nonmatchlist;
}

const vector<alleleValue>& haplotypeCohort::get_alleles_at_site(size_t site_index) const {
  return haplotype_alleles_by_site_index[site_index];
}

vector<size_t> haplotypeCohort::get_active_rows(size_t site, alleleValue a) const {
  if(match_is_rare(site, a)) {
    return get_matches(site, a);
  } else {
    return get_non_matches(site, a);
  }
}

haplotypeCohort::haplotypeCohort(size_t cohort_size, 
            const linearReferenceStructure* reference) : reference(reference) {
  size_t num_sites = reference->number_of_sites();
  number_of_haplotypes = cohort_size;
  haplotype_alleles_by_site_index = vector<vector<alleleValue> >(
              cohort_size,
              vector<alleleValue>());
}

void haplotypeCohort::assign_alleles_at_site(size_t i, 
            vector<alleleValue> alleles_at_site) {
  for(size_t j; j < alleles_at_site.size(); j++) {
    haplotype_alleles_by_site_index[i][j] = alleles_at_site[j];
  }
}