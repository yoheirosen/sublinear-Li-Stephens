#include "reference.hpp"
#include <string.h>
#include <unordered_set>
#include <random>
#include <chrono>
#include <algorithm>

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

size_t linearReferenceStructure::absolute_length() {
  if(final_span_calculated) {
    if(length == 0) {
      length = site_index_to_position.back() + span_lengths.back() + 1;
    }
    return length;
  } else {
    return 0;
  }
}

size_t linearReferenceStructure::absolute_length() const {
  return length;
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

int64_t linearReferenceStructure::add_site(size_t position) {
  if(final_span_calculated) {
    return -3;
  }
  size_t new_index = site_index_to_position.size();
  if(new_index > 0) {
    if(position < site_index_to_position.back()) {
      return -1;
    } else if(position == site_index_to_position.back()) {
      return -2;
    }
    span_lengths.push_back(position - site_index_to_position.back() - 1);
  } else {
    leading_span_length = position - global_offset;
  }
  site_index_to_position.push_back(position);
  position_to_site_index.emplace(position, new_index);
  site_index_to_reference_allele.push_back(unassigned);
  return site_index_to_position.size() - 1;
}

void linearReferenceStructure::add_site(size_t position, 
            alleleValue reference_value) {
  add_site(position);
  site_index_to_reference_allele.back() = reference_value;
}

void linearReferenceStructure::calculate_final_span_length(
            size_t reference_length) {
  size_t previous_position = site_index_to_position.back();
  span_lengths.push_back(reference_length - previous_position - 1);
  final_span_calculated = true;
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
  finalized = true;
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

rowSet haplotypeCohort::get_active_rowSet(size_t site, alleleValue a) const {
  vector<const vector<size_t>* > row_vectors;
  vector<alleleValue> alleles;
  if(match_is_rare(site, a)) {
    row_vectors = {&(haplotype_indices_by_site_and_allele[site][(size_t)a])};
    alleles = {a};
    return rowSet(row_vectors, alleles);
  } else {
    alleleValue b;
    for(size_t i = 0; i < 5; i++) {
      b = (alleleValue)i;
      if(b != a) {
        if(haplotype_indices_by_site_and_allele[site][i].size() != 0) {
          alleles.push_back(b);
          row_vectors.push_back(&(haplotype_indices_by_site_and_allele[site][i]));
        }
      }
    }
    return rowSet(row_vectors, alleles);
  }
}

size_t linearReferenceStructure::pos_ref2global(size_t p) const {
  return p + global_offset;
}

int64_t linearReferenceStructure::pos_global2ref(int64_t p) const {
  return p - global_offset;
}

void linearReferenceStructure::set_allele_at_site(size_t site, alleleValue allele) {
  site_index_to_reference_allele[site] = allele;
}

bool linearReferenceStructure::set_allele_at_pos(size_t p, alleleValue allele) {
  if(is_site(p)) {
    site_index_to_reference_allele[get_site_index(p)] = allele;
    return true;
  } else {
    return false;
  }
}

int linearReferenceStructure::set_allele_at_global_pos(size_t p, alleleValue allele) {
  if(p < global_offset || p >= global_offset + length) {
    return -1;
  } else {
    p = pos_global2ref(p);
    return set_allele_at_pos(p, allele);
  }
}

linearReferenceStructure::linearReferenceStructure(size_t global_offset) : global_offset(global_offset) {
  
}

int haplotypeCohort::add_record(size_t site) {
  if(finalized) {
    return 0;
  } else {
    if(reference->number_of_sites() != 0) {
      if(site == reference->number_of_sites() - 1) {
        haplotype_alleles_by_site_index.push_back(vector<alleleValue>(number_of_haplotypes, unassigned));
        return 1;
      } else {
        return -1;
      }
    } else {
      return -1;
    }
  }
}

int haplotypeCohort::set_sample_allele(size_t site, size_t sample, alleleValue a) {
  if(finalized) {
    return 0;
  } else {
    if(haplotype_alleles_by_site_index[sample][site] == unassigned) {
      haplotype_alleles_by_site_index[sample][site] = a;
      return 1;
    } else {
      return -1;
    }
  }
}

void linearReferenceStructure::set_initial_span(size_t length) {
  leading_span_length = length;
}

void haplotypeCohort::simulate_read_query(
                         const char* ref_seq,
                         double mutation_rate,
                         double recombination_rate,
                         double uncertainty_rate,
                         size_t** return_read_sites,
                         size_t* n_return_read_sites,
                         char** return_read_seq) {
  char* str_to_return = new char[strlen(ref_seq) + 1];
  strcpy(str_to_return, ref_seq);
  return_read_seq = &str_to_return;
  
  // We are making same number of read sites as there are ref sites
  *n_return_read_sites = reference->number_of_sites();
  size_t* read_sites_to_return = new size_t[reference->number_of_sites()];
  
  size_t ref_length = reference->absolute_length();
  
  default_random_engine generator;
  
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  bernoulli_distribution rand_make_uncertain(uncertainty_rate);
  bernoulli_distribution rand_recombine(recombination_rate);
  uniform_int_distribution<size_t> which_haplotype(0, number_of_haplotypes);
  bernoulli_distribution rand_mutate(mutation_rate);
  uniform_int_distribution<size_t> which_allele(0, 4);
  uniform_int_distribution<size_t> which_position(0, ref_length - 1);
  
  // Deal with sites ***********************************************************
  
  size_t ref_sites = reference->number_of_sites();
  size_t read_sites_left = ref_sites;
  unordered_set<size_t> blacklist;
  
  vector<size_t> read_sites;
  
  // decide which sites are shared
  for(size_t i = 0; i < ref_sites; i++) {
    blacklist.emplace(reference->get_position(i));
    if(rand_make_uncertain(generator)) {
      read_sites.push_back(reference->get_position(i));
      --read_sites_left;
    }
  }
  
  // decide which non-ref sites are read sites
  while(read_sites_left > 0) {
    size_t candidate = which_position(generator);
    if(blacklist.count(candidate) == 0) {
      read_sites.push_back(candidate);
      blacklist.emplace(candidate);
      --read_sites_left;
    }
  }
  
  // need to return ordered size_t C-style array
  sort(read_sites.begin(), read_sites.end());
  for(size_t i = 0; i < read_sites.size(); i++) {
    read_sites_to_return[i] = read_sites[i];
  }
  
  // Deal with alleles *********************************************************
  
  size_t pos;
  for(pos = 0; pos < reference->get_position(0); pos++) {
    if(blacklist.count(pos) == 0) {
      if(rand_mutate(generator)) {
        size_t mutate_to = which_allele(generator);
        str_to_return[pos] = allele_to_char((alleleValue)mutate_to);
      }
    }
  }

  // will piece together haplotype using Li-Stephens approximation
  size_t h_index = which_haplotype(generator);
  for(size_t site = 0; site < ref_sites; site++) {
    str_to_return[reference->get_position(site)] = allele_to_char(allele_at(site, h_index));
    if(rand_recombine(generator)) {
      h_index = which_haplotype(generator);
    }
    for(size_t pos = reference->get_position(site) + 1; pos < reference->get_position(site + 1); pos++) {
      if(rand_mutate(generator)) {
        size_t mutate_to = which_allele(generator);
        str_to_return[pos] = allele_to_char((alleleValue)mutate_to);
      }
      if(rand_recombine(generator)) {
        h_index = which_haplotype(generator);
      }
    }
  }  
}