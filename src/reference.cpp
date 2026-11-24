#include "reference.hpp"
#include <string.h>
#include <unordered_set>
#include <random>
#include <chrono>
#include <algorithm>
#include <list>
#include <htslib/vcf.h>
#include <iostream>

using namespace std;

siteIndex::siteIndex(size_t global_offset) : global_offset(global_offset) {
}

siteIndex::siteIndex(const vector<string>& haplotypes) {
  size_t number_of_haplotypes = haplotypes.size();
  size_t length_of_haplotypes = haplotypes[0].size();
  for(size_t i = 0; i < length_of_haplotypes; i++) {
    vector<alleleValue> alleles_seen;
    for(size_t j = 0; j < number_of_haplotypes; j++) {
      alleleValue allele_at_j = allele::from_char(haplotypes[j][i]);
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
      add_site(i);
    }
  }
  calculate_final_span_length(length_of_haplotypes);
}

siteIndex::siteIndex(const vector<size_t>& positions, size_t length) {
  for(size_t i = 0; i < positions.size(); i++) {
    add_site(positions[i]);
  }
  calculate_final_span_length(length);
}

siteIndex::siteIndex(const vector<size_t>& positions, size_t length, size_t global_offset) :
  site_index_to_position(positions), length(length), global_offset(global_offset) {
  calculate_final_span_length(length);
}

siteIndex::~siteIndex() {
  
}

size_t haplotypeCohort::get_n_sites() const {
  return allele_counts_by_site_index.size();
}

bool siteIndex::is_site(size_t actual_position) const {
  return (position_to_site_index.count(actual_position) == 1);
}

size_t siteIndex::get_site_index(size_t actual_position) const {
  return position_to_site_index.at(actual_position);
}

size_t siteIndex::get_position(size_t site_index) const {
  return site_index_to_position[site_index];
}

bool siteIndex::has_span_before(size_t site_index) const {
  if(site_index == 0) {
    return (leading_span_length != 0);
  } else {
    return (span_lengths[site_index - 1] != 0);
  }
}

bool siteIndex::has_span_after(size_t site_index) const {
  return (span_lengths[site_index] != 0);
}

size_t siteIndex::span_length_before(size_t site_index) const {
  if(site_index == 0) {
    return leading_span_length;
  } else {
    return span_lengths[site_index - 1];
  }
}

size_t siteIndex::span_length_after(size_t site_index) const {
  return span_lengths[site_index];
}

size_t siteIndex::number_of_sites() const {
  return site_index_to_position.size();
}

size_t siteIndex::length_in_bp() const {
  return length;
}

int64_t siteIndex::add_site(size_t position) {
  size_t new_index = site_index_to_position.size();
  if(new_index > 0) {
    if(position < site_index_to_position.back()) {
      throw runtime_error("constructed siteIndex out of order");
    } else if(position == site_index_to_position.back()) {
      throw runtime_error("double-wrote siteIndex site");
    }
    span_lengths.push_back(position - site_index_to_position.back() - 1);
  } else {
    leading_span_length = position - global_offset;
  }
  site_index_to_position.push_back(position);
  position_to_site_index.emplace(position, new_index);
  return site_index_to_position.size() - 1;
}

void siteIndex::calculate_final_span_length(size_t reference_length) {
  if(site_index_to_position.size() > 0) {
    size_t previous_position = site_index_to_position.back();
    span_lengths.push_back(global_offset + reference_length - previous_position - 1);
    length = site_index_to_position.back() + span_lengths.back() + 1 - global_offset;
  }
}

size_t siteIndex::find_site_above(size_t position) const {
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

size_t siteIndex::find_site_below(size_t position) const {
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
  return alleles_by_haplotype_and_site[haplotype_index][site_index];
}

const vector<alleleValue>& haplotypeCohort::get_haplotype(size_t idx) const {
  return alleles_by_haplotype_and_site[idx];
}

size_t haplotypeCohort::number_matching(size_t site_index, alleleValue a) const {
  return allele_counts_by_site_index[site_index][(int)a];
}

size_t haplotypeCohort::number_not_matching(size_t site_index, alleleValue a) const {
  return number_of_haplotypes - number_matching(site_index, a);
}

size_t haplotypeCohort::number_active(site_idx_t site_index, alleleValue a) const {
  return match_is_rare(site_index, a) ? number_matching(site_index, a) : number_not_matching(site_index, a);
}

void haplotypeCohort::populate_allele_counts() {
  int num_sites = reference->number_of_sites();
  allele_counts_by_site_index = vector<vector<size_t> >(num_sites,
              vector<size_t>(5, 0));
  haplotype_indices_by_site_and_allele = vector<vector<vector<size_t> > >(
              num_sites, vector<vector<size_t> >(5, vector<size_t>()));
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    for(int j = 0; j < num_sites; j++) {
      int allele_rank = (int)alleles_by_haplotype_and_site[i][j];
      allele_counts_by_site_index[j][allele_rank]++;
      haplotype_indices_by_site_and_allele[j][allele_rank].push_back(i);
    }
  }

  active_rowSets_by_site_and_allele = vector<vector<rowSet> >(num_sites, vector<rowSet>(5));
  for(size_t i = 0; i < num_sites; i++) {
    for(size_t a = 0; a < 5; a++) {
      if(haplotype_indices_by_site_and_allele[i][a].size() > number_of_haplotypes / 2) {
        haplotype_indices_by_site_and_allele[i][a].clear();
      }
      active_rowSets_by_site_and_allele[i][a] = build_active_rowSet(i, (alleleValue)a);
    }
  }
  finalized = true;
}

size_t haplotypeCohort::get_n_haplotypes() const {
  return number_of_haplotypes;
}

// size_t haplotypeCohort::get_data_size() {
//   size_t total_size = get_n_sites() * (5 * (sizeof(size_t) + (5 * sizeof(const vector<size_t>*)))) + sizeof(siteIndex*) + sizeof(size_t) + sizeof(bool);
//   for(size_t i = 0; i < get_n_sites; i++) {
//     for(size_t a = 0; a < 5; a++) {
//       if(allele_counts_by_site_index[i][a] <= number_of_haplotypes / 2) {
//         total_size += sizeof(size_t) * allele_counts_by_site_index[i][a];
//       }
//     }
//   }
//   return total_size;
// }
// 
// size_t siteIndex::get_data_size() {
//   return ((2 * number_of_sites()) + 3) * sizeof(size_t);
// }


const vector<size_t>& haplotypeCohort::get_matches(size_t site_index, alleleValue a) const {
  size_t allele_rank = (size_t)a;
  return haplotype_indices_by_site_and_allele[site_index][allele_rank];
}

haplotypeCohort::haplotypeCohort(const vector<string>& haplotypes, 
            siteIndex* reference) : reference(reference) {
  number_of_haplotypes = haplotypes.size();
  size_t num_sites = reference->number_of_sites();
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    alleles_by_haplotype_and_site.push_back(vector<alleleValue>(num_sites));
    for(size_t j = 0; j < num_sites; j++) {
      size_t pos = reference->get_position(j);
      alleles_by_haplotype_and_site[i][j] = allele::from_char(haplotypes[i][pos],
                  unassigned);
    }
  }
  populate_allele_counts();
}

haplotypeCohort::haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
          siteIndex* reference) : reference(reference),
          alleles_by_haplotype_and_site(haplotypes) {
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

// const vector<alleleValue>& haplotypeCohort::allele_vector_at_site(size_t site_index) const {

// }

vector<size_t> haplotypeCohort::get_active_rows(size_t site, alleleValue a) const {
  if(match_is_rare(site, a)) {
    return get_matches(site, a);
  } else {
    return get_non_matches(site, a);
  }
}

haplotypeCohort::haplotypeCohort(size_t cohort_size, 
            siteIndex* reference) : reference(reference) {
  size_t num_sites = reference->number_of_sites();
  number_of_haplotypes = cohort_size;
  alleles_by_haplotype_and_site = vector<vector<alleleValue> >(
              cohort_size,
              vector<alleleValue>(num_sites, unassigned));
}

void haplotypeCohort::assign_alleles_at_site(size_t i, 
            vector<alleleValue> alleles_at_site) {
  for(size_t j = 0; j < alleles_at_site.size(); j++) {
    alleles_by_haplotype_and_site[i][j] = alleles_at_site[j];
  }
}

const rowSet& haplotypeCohort::get_active_rowSet(size_t site, alleleValue a) const {
  return active_rowSets_by_site_and_allele[site][(size_t)a];
}

rowSet haplotypeCohort::build_active_rowSet(size_t site, alleleValue a) const {
  if(allele_counts_by_site_index[site][(size_t)a] == 0 || allele_counts_by_site_index[site][(size_t)a] == number_of_haplotypes) {
    return rowSet();
  } else {
    vector<const vector<size_t>* > row_vectors;
    if(match_is_rare(site, a)) {
      row_vectors = {&(haplotype_indices_by_site_and_allele[site][(size_t)a])};
      return rowSet(row_vectors);
    } else {
      alleleValue b;
      for(size_t i = 0; i < 5; i++) {
        b = (alleleValue)i;
        if(b != a) {
          if(haplotype_indices_by_site_and_allele[site][i].size() != 0) {
            row_vectors.push_back(&(haplotype_indices_by_site_and_allele[site][i]));
          }
        }
      }
      return rowSet(row_vectors);
    }
  }
}

size_t siteIndex::pos_ref2global(size_t p) const {
  return p + global_offset;
}

int64_t siteIndex::pos_global2ref(int64_t p) const {
  return p - global_offset;
}

size_t siteIndex::start_position() const {
  return global_offset;
}

size_t siteIndex::end_position() const {
  return start_position() + length_in_bp() - 1;
}

void haplotypeCohort::add_record() {
  if(finalized) {
    throw runtime_error("attempted to add record to locked haplotype cohort");
  } else {
    if(reference->number_of_sites() != 0 && get_n_sites() < reference->number_of_sites()) {
      for(size_t i = 0; i < number_of_haplotypes; i++) {
        alleles_by_haplotype_and_site[i].push_back(unassigned);
      }
      return;
    } else {
      throw runtime_error("attempted to add more records than number of sites in reference");
    }
  }
}

void haplotypeCohort::set_sample_allele(size_t site, size_t sample, alleleValue a) {
  if(finalized) {
    throw runtime_error("attempted to modify locked haplotype cohort");
  } else {
    if(alleles_by_haplotype_and_site[sample][site] == unassigned) {
      alleles_by_haplotype_and_site[sample][site] = a;
      return;
    } else {
      throw runtime_error("attempted to double-write haplotype allele");
    }
  }
}

void siteIndex::set_initial_span(size_t length) {
  leading_span_length = length;
}

alleleValue haplotypeCohort::get_dominant_allele(size_t site) const {
  size_t max_count = 0;
  alleleValue dominant_allele = unassigned;
  for(size_t i = 0; i < allele_counts_by_site_index[site].size(); i++) {
    if(allele_counts_by_site_index[site][i] > max_count) {
      max_count = allele_counts_by_site_index[site][i];
      dominant_allele = (alleleValue)i;
    }
  }
  return dominant_allele;
}

size_t haplotypeCohort::get_total_information(size_t site) const {
  size_t to_return = 0;
  for(size_t a = 0; a < 5; a++) {
    to_return += haplotype_indices_by_site_and_allele[site][(size_t)a].size();
  }
  return to_return;
}

size_t haplotypeCohort::get_information_content(size_t site, alleleValue a) const {
  if(match_is_rare(site, a)) {
    return number_matching(site, a);
  } else {
    return number_of_haplotypes - number_matching(site, a);
  }
}

size_t haplotypeCohort::sum_information_content(const vector<alleleValue>& query, size_t start_site) const {
  size_t sum = 0;
  for(size_t i = 0; i < query.size(); i++) {
    sum += get_information_content(i + start_site, query[i]);
  }
  return sum;
}

void haplotypeCohort::set_column(const vector<alleleValue>& values) {
  set_column(values, get_n_sites() - 1);
}

void haplotypeCohort::set_column(const vector<alleleValue>& values, size_t site) {
  for(size_t i = 0; i < get_n_haplotypes(); i++) {
    alleles_by_haplotype_and_site[i][site] = values[i];
  }
}

siteIndex* haplotypeCohort::get_reference() const {
  return reference;
}

haplotypeCohort* haplotypeCohort::remove_sites_below_frequency(double frequency) const {
  size_t biggest_major = get_n_sites() * (1 - frequency);
  vector<size_t> remaining_sites; 
  for(size_t i = 0; i < get_n_sites(); i++) {
    bool passes = true;
    for(size_t j = 0; j < 5; j++) {
      if(allele_counts_by_site_index[i][j] > biggest_major) {
        passes = false;
      }
    }
    if(passes) {
      remaining_sites.push_back(i);
    }
  }
  vector<size_t> remaining_site_positions(remaining_sites.size());
  for(size_t i = 0; i < remaining_sites.size(); i++) {
    remaining_site_positions[i] = reference->get_position(remaining_sites[i]);
  }
  siteIndex* new_reference = 
            new siteIndex(remaining_site_positions,
                                         reference->length_in_bp(),
                                         reference->start_position() +
                                                  reference->length_in_bp());
  haplotypeCohort* to_return = 
            new haplotypeCohort(number_of_haplotypes, new_reference);
  for(size_t i = 0; i < remaining_sites.size(); i++) {
    to_return->add_record();
    to_return->set_column(alleles_by_haplotype_and_site[remaining_sites[i]]);
  }
  to_return->populate_allele_counts();
  return to_return;
}

vector<size_t> haplotypeCohort::rand_haplos(size_t N) const {
  return haploRandom::n_unique_uints(N, get_n_haplotypes());
}

vector<size_t> siteIndex::rand_sites(size_t N) const {
  return haploRandom::n_unique_uints(N, number_of_sites());
}

vector<size_t> siteIndex::rand_site_positions(size_t N) const {
  return haploRandom::n_unique_uints(N, number_of_sites(), &site_index_to_position);
}

size_t siteIndex::rand_interval_start(size_t len) const {
  if(len > length_in_bp()) {
    throw runtime_error("attempted to generate random subinterval exceeding size of reference");
  }
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> unit_uniform(0.0, 1.0);
  size_t supremum = length_in_bp() - len;
  size_t draw = (size_t)(unit_uniform(generator) * supremum);
  if(draw == supremum) { draw = supremum - 1; }
  return draw + global_offset;
}

size_t siteIndex::rand_site_interval_start(size_t N) const {
  if(N > number_of_sites()) {
    throw runtime_error("attempted to generate random subinterval exceeding size of reference");
  }
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> unit_uniform(0.0, 1.0);
  size_t supremum = number_of_sites() - N;
  size_t draw = (size_t)(unit_uniform(generator) * supremum);
  if(draw == supremum) { draw = supremum - 1; }
  return draw;
}

size_t haplotypeCohort::rand_haplo_idx() const {
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> unit_uniform(0.0, 1.0);
  size_t supremum = get_n_haplotypes();
  size_t draw = (size_t)(unit_uniform(generator) * supremum);
  if(draw == supremum) { draw = supremum - 1; }
  return draw;
}

size_t haplotypeCohort::rand_haplo_idx(size_t current) const {
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> unit_uniform(0.0, 1.0);
  size_t supremum = get_n_haplotypes() - 1;
  size_t draw = (size_t)(unit_uniform(generator) * supremum);
  if(draw == supremum) { draw = supremum - 1; }
  if(draw >= current) { draw++; }
  return draw;
}

vector<alleleValue> haplotypeCohort::rand_LS_haplo(double log_recomb_probability, double log_mutation_probability) const {
  return rand_LS_haplo(log_recomb_probability, log_mutation_probability, 0, get_n_sites() - 1);
}

vector<alleleValue> haplotypeCohort::rand_LS_haplo(double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const {
  if(end_site > reference->number_of_sites() || start_site > end_site) {
    throw runtime_error("invalid bounds for random haplotype generation");
  }
  vector<alleleValue> to_return(end_site - start_site + 1);
  size_t h_idx = rand_haplo_idx();
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  for(size_t i = 0; i <= end_site - start_site - 1; i++) {
    binomial_distribution<size_t> recombiner(get_reference()->span_length_after(i) + 1, exp(log_recomb_probability));
    size_t n_recombinations = recombiner(generator);
    to_return[i] = haploRandom::mutate(allele_at(i + start_site, h_idx), log_mutation_probability);
    for(size_t j = 0; j < n_recombinations; j++) {
      h_idx = rand_haplo_idx(h_idx);
    }
  }
  to_return.back() = haploRandom::mutate(allele_at(end_site - 1, h_idx), log_mutation_probability);
  return to_return;
}

vector<alleleValue> haplotypeCohort::rand_desc_haplo(size_t generations, double log_recomb_probability, double log_mutation_probability) const {
  return rand_desc_haplo(generations, log_recomb_probability, log_mutation_probability, 0, get_n_sites() - 1);
}

vector<alleleValue> haplotypeCohort::rand_desc_haplo(size_t generations, double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const {
  vector<vector<alleleValue> > old_gen;
  vector<vector<alleleValue> > new_gen;
  
  size_t gen_size = (size_t)pow(2, generations);
  new_gen.clear();
  for(size_t j = 0; j < gen_size; j++) {
    new_gen.push_back(rand_LS_haplo(log_recomb_probability, log_mutation_probability, start_site, end_site));
  }
  
  for(size_t i = 1; i < generations; i++) {
    gen_size = (size_t)pow(2, generations - i);
    old_gen = new_gen;
    new_gen.clear();
    for(size_t j = 0; j < gen_size; j++) {
      new_gen.push_back(get_reference()->make_child(old_gen[2*j], old_gen[2*j + 1], log_recomb_probability, log_mutation_probability, start_site, end_site));
    }
  }
  
  return(get_reference()->make_child(old_gen[0], old_gen[1], log_recomb_probability, log_mutation_probability, start_site, end_site));
}

vector<alleleValue> siteIndex::make_child(const vector<alleleValue>& parent_0, const vector<alleleValue>& parent_1, double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const {
  vector<alleleValue> to_return(end_site - start_site + 1);
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  bernoulli_distribution hap_chooser(0.5);
  bool h_idx = hap_chooser(generator);
  
  for(size_t i = 0; i <= end_site - start_site - 1; i++) {
    binomial_distribution<size_t> recombiner(span_length_after(i + start_site) + 1, exp(log_recomb_probability));
    size_t n_recombinations = recombiner(generator);
    to_return[i] = h_idx ? parent_0[i] : parent_1[i];
    to_return[i] = haploRandom::mutate(to_return[i], log_mutation_probability);
    h_idx = (bool)(((int)h_idx + n_recombinations) % 2);
  }
  to_return.back() = h_idx ? parent_0.back() : parent_1.back();
  to_return.back() = haploRandom::mutate(to_return.back(), log_mutation_probability);
  return to_return;
}

vector<size_t> haploRandom::n_unique_uints(size_t N, size_t supremum) {
  return haploRandom::n_unique_uints(N, supremum, nullptr);
}

vector<size_t> haploRandom::n_unique_uints(size_t N, size_t supremum, const vector<size_t>* blacklist) {
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  uniform_real_distribution<double> unit_uniform(0.0, 1.0);
  vector<size_t> draws(N);
  for(size_t i = 0; i < N; i++) {
    size_t draw = (size_t)(unit_uniform(generator) * supremum);
    if(draw == supremum) { draw = supremum - 1; } // mathematically valid assuming true real numbers; good enough here
    draws[i] = draw;
    supremum--;
  }
  
  list<size_t> actual_values;
  
  if(blacklist != nullptr) {
    for(size_t i = 0; i < blacklist->size(); i++) {
      actual_values.push_back(blacklist->at(i));
    }
    actual_values.sort();
  }
  
  for(size_t i = 0; i < draws.size(); i++) {
    auto it = actual_values.begin();
    size_t j;
    while(it != actual_values.end()) {
      if(*it <= draws[i]) {
        draws[i]++;
        ++it;
      } else {
        break;
      }
    }
    actual_values.insert(it, draws[i]);
  }
  
  if(blacklist != nullptr) {
    for(size_t i = 0; i < blacklist->size(); i++) {
      actual_values.remove(blacklist->at(i));
    }
  }
  
  vector<size_t> to_return;
  
  auto it = actual_values.begin();
  do {
    to_return.push_back(*it);
    ++it;
  } while (it != actual_values.end());
  
  return to_return;
}

alleleValue haploRandom::mutate(alleleValue a, double log_mutation_probability) {
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  bernoulli_distribution mutate(exp(log_mutation_probability));
  if(mutate(generator)) {
    uniform_real_distribution<double> unit_uniform(0.0, 1.0);
    size_t supremum = 4;
    size_t draw = (size_t)(unit_uniform(generator) * supremum);
    if(draw == supremum) { draw = supremum - 1; }
    if(draw >= (size_t)a) { ++draw; }
    return (alleleValue)draw;
  }
  return a;
}

haplotypeCohort* haplotypeCohort::remove_rare_sites(double max_rarity) const {
  size_t maj_all_fq_limit = (size_t)(get_n_haplotypes() * (1 - max_rarity));
  vector<size_t> sites_not_dropped; 
  for(size_t i = 0; i < get_n_sites(); i++) {
    bool passes = true;
    for(size_t j = 0; j < 5; j++) {
      if(allele_counts_by_site_index[i][j] > maj_all_fq_limit) {
        passes = false;
      }
    }
    if(passes) {
      sites_not_dropped.push_back(i);
    }
  }
  vector<size_t> remaining_site_positions(sites_not_dropped.size());
  for(size_t i = 0; i < sites_not_dropped.size(); i++) {
    remaining_site_positions[i] = reference->get_position(sites_not_dropped[i]);
  }
  siteIndex* new_ref = new siteIndex(
                       remaining_site_positions,
                       reference->length_in_bp(),
                       reference->start_position() + reference->length_in_bp());
  haplotypeCohort* to_return = 
            new haplotypeCohort(number_of_haplotypes, new_ref);
  for(size_t i = 0; i < sites_not_dropped.size(); i++) {
    vector<alleleValue> this_column(number_of_haplotypes, unassigned);
    for(size_t j = 0; j < number_of_haplotypes; j++) {
      this_column[j] = alleles_by_haplotype_and_site[j][sites_not_dropped[i]];
    }
    to_return->set_column(this_column, i);
  }
  to_return->populate_allele_counts();
  return to_return;
}

void haplotypeCohort::remove_homogeneous_sites() {
  vector<size_t> sites_to_keep;
  for(size_t i = 0; i < get_n_sites(); i++) {
    size_t n_alleles_seen = 0;
    for(size_t a = 0; a < 5; a++) {
      if(number_matching(i, (alleleValue)a) != 0) {
        ++n_alleles_seen;
      }
    }
    if(n_alleles_seen > 1) {
      sites_to_keep.push_back(i);
    }
  }
  
  reference->keep_subset_of_sites(sites_to_keep);
  
  if(sites_to_keep.size() != 0) {
    for(size_t i = 0; i < sites_to_keep.size(); i++) {
      size_t old_i = sites_to_keep[i];
      for(size_t j = 0; j < number_of_haplotypes; j++) {
        alleles_by_haplotype_and_site[j][i] = alleles_by_haplotype_and_site[j][old_i];
      }
      haplotype_indices_by_site_and_allele[i] = haplotype_indices_by_site_and_allele[old_i];
      allele_counts_by_site_index[i] = allele_counts_by_site_index[old_i];
    }
    for(size_t j = 0; j < number_of_haplotypes; j++) {
      alleles_by_haplotype_and_site[j].resize(sites_to_keep.size());
    }
    haplotype_indices_by_site_and_allele.resize(sites_to_keep.size());
    allele_counts_by_site_index.resize(sites_to_keep.size());
    
    active_rowSets_by_site_and_allele = vector<vector<rowSet> >(sites_to_keep.size(), vector<rowSet>(5));
    for(size_t i = 0; i < sites_to_keep.size(); i++) {
      for(size_t a = 0; a < 5; a++) {
        active_rowSets_by_site_and_allele[i][a] = build_active_rowSet(i, (alleleValue)a);
      }
    }
  } else {
    alleles_by_haplotype_and_site.clear();
    haplotype_indices_by_site_and_allele.clear();
    active_rowSets_by_site_and_allele.clear();
    allele_counts_by_site_index.clear();
  }
}

void siteIndex::keep_subset_of_sites(const vector<size_t>& sites_to_keep) {
  if(sites_to_keep.size() > 0) {
    vector<size_t> new_span_lengths(sites_to_keep.size());
    size_t new_leading_span;
    vector<size_t> new_positions(sites_to_keep.size());
    for(size_t i = 0; i < sites_to_keep.size(); i++) {
      new_positions[i] = get_position(sites_to_keep[i]);
    }
    new_leading_span = new_positions[0] - global_offset;
    for(size_t i = 0; i < sites_to_keep.size() - 1; i++) {
      new_span_lengths[i] = new_positions[i + 1] - new_positions[i] - 1;
    }
    new_span_lengths[sites_to_keep.size() - 1] = global_offset + length - new_positions[sites_to_keep.size() - 1] - 1;
    leading_span_length = new_leading_span;
    span_lengths = new_span_lengths;
    site_index_to_position = new_positions;
    {
      unordered_map<size_t, size_t> temp;
      std::swap(position_to_site_index, temp);
    }
    for(size_t i = 0; i < new_positions.size(); i++) {
      position_to_site_index.emplace(new_positions[i], i);
    }
  } else {
    leading_span_length = length;
    position_to_site_index.clear();
    site_index_to_position.clear();
    span_lengths.clear();
  }
}

haplotypeCohort* haplotypeCohort::subset(size_t start_site, size_t end_site, size_t n_to_keep) const {
  vector<size_t> ids_to_keep = rand_haplos(n_to_keep);
  return subset(start_site, end_site, ids_to_keep);
}

haplotypeCohort* haplotypeCohort::subset(size_t start_site, size_t end_site, const vector<size_t>& ids_to_keep) const {
  size_t n_sites = end_site - start_site + 1;
  vector<size_t> remaining_site_positions(n_sites);
  for(size_t i = 0; i < n_sites; i++) {
    remaining_site_positions[i] = reference->get_position(i + start_site);
  }
  siteIndex* new_ref = new siteIndex(remaining_site_positions, reference->length_in_bp(), reference->start_position());
  
  vector<vector<alleleValue> > kept_haplotypes(ids_to_keep.size(), vector<alleleValue>(n_sites));
  
  for(size_t i = start_site; i < start_site + n_sites; i++) {
    for(size_t j_it = 0; j_it < ids_to_keep.size(); j_it++) {
      size_t j = ids_to_keep[j_it];
      kept_haplotypes[j_it][i - start_site] = allele_at(i, j);
    }
  }
  
  haplotypeCohort* to_return = new haplotypeCohort(kept_haplotypes, new_ref);
  to_return->populate_allele_counts();
  to_return->remove_homogeneous_sites();
  return to_return;
}

haplotypeCohort* build_cohort(const string& vcf_path) {
  vcfFile* cohort_vcf = vcf_open(vcf_path.c_str(), "r");
  bcf_hdr_t* cohort_hdr = bcf_hdr_read(cohort_vcf);
  bcf1_t* record = bcf_init1();
  
  size_t number_of_haplotypes = bcf_hdr_nsamples(cohort_hdr) * 2;
  
  siteIndex* reference = new siteIndex(0); 
  haplotypeCohort* cohort = new haplotypeCohort(number_of_haplotypes, reference);
  bool built_initial_span = false;
  bool at_start = true;
  size_t last_site_position = 0;
  while(bcf_read(cohort_vcf, cohort_hdr, record) == 0) {
    size_t site_position = record->pos;
    if(bcf_is_snp(record) == 1) {
      if(site_position <= last_site_position && last_site_position != 0 && !at_start) {
        // stderr << "skipping overlapping position " << site_position << " which overlaps " << last_site_position << endl;
      } else {
        last_site_position = site_position;
        
        bcf_unpack(record, BCF_UN_ALL);
        if(!built_initial_span) {
          reference->set_initial_span(site_position);
          built_initial_span = true;
        }
        
        int64_t site_index = reference->add_site(site_position);
        cohort->add_record();
        
        int32_t *gt_arr = NULL, ngt_arr = 0;
        int ngt = bcf_get_genotypes(cohort_hdr, record, &gt_arr, &ngt_arr);
        if(site_index >= 0) {
          for(size_t i = 0; i < number_of_haplotypes; i++) {
            int allele_index = bcf_gt_allele(gt_arr[i]);
            char allele_value = record->d.allele[allele_index][0];
            cohort->set_sample_allele(site_index, i, allele::from_char(allele_value));
          }
        }
        free(gt_arr);
      }
      at_start = false;
    }
  }
  
  size_t ref_end = last_site_position;
  
  // cerr << "loaded vcf " << vcf_path << endl;
  reference->calculate_final_span_length(last_site_position);
  cohort->populate_allele_counts();
  // cerr << "built haplotypecohort object" << endl;
  
  bcf_hdr_destroy(cohort_hdr);
  bcf_destroy(record);
  vcf_close(cohort_vcf);
  
  return cohort;
}

void siteIndex::serialize_human(std::ostream& indexout) const {
  indexout << global_offset << "\t" << length << "\t" << site_index_to_position.size() << endl;
  for(size_t i = 0; i < site_index_to_position.size(); i++) {
    indexout << site_index_to_position[i] << "\t";
  }
  indexout << endl;
}

void haplotypeCohort::serialize_human(std::ostream& cohortout) const {
  reference->serialize_human(cohortout);
  cohortout << number_of_haplotypes << endl;
  for(size_t i = 0; i < get_n_sites(); i++) {
    for(size_t a = 0; a < 5; a++) {
      cohortout << allele_counts_by_site_index[i][a] << "\t";
    }
    cohortout << endl;
    for(size_t a = 0; a < 5; a++) {
      for(size_t j = 0; j < haplotype_indices_by_site_and_allele[i][a].size(); j++) {
        cohortout << haplotype_indices_by_site_and_allele[i][a][j] << "\t";
      }
    }
    cohortout << endl; 
  }
}

siteIndex::siteIndex(std::istream& indexin) {
  indexin >> global_offset;
  indexin >> length;
  size_t n_sites;
  indexin >> n_sites;
  site_index_to_position = vector<size_t>(n_sites);
  span_lengths = vector<size_t>(n_sites);
  for(size_t i = 0; i < n_sites; i++) {
    indexin >> site_index_to_position[i];
  }
  leading_span_length = site_index_to_position[0] - global_offset;
  for(size_t i = 1; i < n_sites - 1; i++) {
    span_lengths[i - 1] = site_index_to_position[i] - site_index_to_position[i - 1] - 1;
  }
  span_lengths[n_sites - 1] = global_offset + length - site_index_to_position[n_sites - 1] - 1;
  for(size_t i = 0; i < n_sites; i++) {
    position_to_site_index.emplace(site_index_to_position[i], i);
  }
}

haplotypeCohort::haplotypeCohort(std::istream& cohortin, siteIndex* reference) : reference(reference) {
  cohortin >> number_of_haplotypes;
  haplotype_indices_by_site_and_allele = vector<vector<vector<haplo_id_t> > >(reference->number_of_sites(),
                                                vector<vector<haplo_id_t> >(5));
  active_rowSets_by_site_and_allele = vector<vector<rowSet> >(reference->number_of_sites(),
                                             vector<rowSet>(5));
  allele_counts_by_site_index = vector<vector<size_t> >(reference->number_of_sites(),
                                       vector<size_t>(5));
  for(size_t i = 0; i < reference->number_of_sites(); i++) {
    for(size_t a = 0; a < 5; a++) {
      cohortin >> allele_counts_by_site_index[i][a];
      if(allele_counts_by_site_index[i][a] <= number_of_haplotypes / 2) {
        haplotype_indices_by_site_and_allele[i][a] = vector<haplo_id_t>(allele_counts_by_site_index[i][a]);
      } else {
        haplotype_indices_by_site_and_allele[i][a] = vector<haplo_id_t>(0);
      }
    }
    for(size_t a = 0; a < 5; a++) {
      for(size_t j = 0; j < haplotype_indices_by_site_and_allele[i][a].size(); j++) {
        cohortin >> haplotype_indices_by_site_and_allele[i][a][j];
      }
    }
    for(size_t a = 0; a < 5; a++) {
      active_rowSets_by_site_and_allele[i][a] = build_active_rowSet(i, (alleleValue)a);
    }
  }
  finalized = true;
}

void haplotypeCohort::uncompress() {
  size_t n_sites = get_n_sites();
  if(alleles_by_haplotype_and_site.size() == 0) {
    vector<alleleValue> prototype = vector<alleleValue>(n_sites, unassigned);
    for(size_t i = 0; i < n_sites; i++) {
      for(size_t a = 0; a < 5; a++) {
        if(!match_is_rare(i, (alleleValue)a)) {
          prototype[i] = (alleleValue)a;
        }
      }
    }
    alleles_by_haplotype_and_site = vector<vector<alleleValue>>(number_of_haplotypes, prototype);
    for(size_t i = 0; i < n_sites; i++) {
      for(size_t a = 0; a < 5; a++) {
        if((alleleValue)a != prototype[i]) {
          alleleValue allele = (alleleValue)a;
          for(size_t j_it = 0; j_it < haplotype_indices_by_site_and_allele[i][a].size(); j_it++) {
            alleles_by_haplotype_and_site[haplotype_indices_by_site_and_allele[i][a][j_it]][i] = allele;
          }
        }
      }
    }
  }
}

void haplotypeCohort::compress() {
  alleles_by_haplotype_and_site.clear();
  alleles_by_haplotype_and_site.shrink_to_fit();
  for(size_t i = 0; i < get_n_sites(); i++) {
    for(size_t a = 0; a < 5; a++) {
      if(!match_is_rare(i, (alleleValue)a)) {
        haplotype_indices_by_site_and_allele[i][a].clear();
        haplotype_indices_by_site_and_allele[i][a].shrink_to_fit();
      }
    }
  }
}