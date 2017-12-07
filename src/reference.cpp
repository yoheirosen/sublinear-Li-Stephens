#include "reference.hpp"
#include <string.h>
#include <unordered_set>
#include <random>
#include <chrono>
#include <algorithm>
#include <list>

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

size_t siteIndex::absolute_length() const {
  return length;
}

int64_t siteIndex::add_site(size_t position) {
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
  return site_index_to_position.size() - 1;
}

void siteIndex::calculate_final_span_length(
            size_t reference_length) {
  size_t previous_position = site_index_to_position.back();
  span_lengths.push_back(reference_length - previous_position - 1);
  length = site_index_to_position.back() + span_lengths.back() + 1;
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

size_t haplotypeCohort::get_n_haplotypes() const {
  return number_of_haplotypes;
}

const vector<size_t>& haplotypeCohort::get_matches(size_t site_index, alleleValue a) const {
  size_t allele_rank = (size_t)a;
  return haplotype_indices_by_site_and_allele[site_index][allele_rank];
}

haplotypeCohort::haplotypeCohort(const vector<string>& haplotypes, 
            const siteIndex* reference) : reference(reference) {
  number_of_haplotypes = haplotypes.size();
  size_t num_sites = reference->number_of_sites();
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    haplotype_alleles_by_site_index.push_back(vector<alleleValue>(num_sites));
    for(size_t j = 0; j < num_sites; j++) {
      size_t pos = reference->get_position(j);
      haplotype_alleles_by_site_index[i][j] = allele::from_char(haplotypes[i][pos],
                  unassigned);
    }
  }
  populate_allele_counts();
}

haplotypeCohort::haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
          const siteIndex* reference) : reference(reference),
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

const vector<alleleValue>& haplotypeCohort::allele_vector_at_site(size_t site_index) const {
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
            const siteIndex* reference) : reference(reference) {
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

size_t siteIndex::pos_ref2global(size_t p) const {
  return p + global_offset;
}

int64_t siteIndex::pos_global2ref(int64_t p) const {
  return p - global_offset;
}

size_t siteIndex::start_position() const {
  return global_offset;
}


int haplotypeCohort::add_record(size_t site) {
  if(finalized) {
    return 0;
  } else {
    if(reference->number_of_sites() != 0) {
      if(site == reference->number_of_sites() - 1) {
        for(size_t i = 0; i < number_of_haplotypes; i++) {
          haplotype_alleles_by_site_index[i].push_back(unassigned);
        }
        // haplotype_alleles_by_site_index.push_back(vector<alleleValue>(number_of_haplotypes, unassigned));
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

size_t haplotypeCohort::get_MAC(size_t site) const {
  return number_of_haplotypes - number_matching(site, get_dominant_allele(site));
}

size_t haplotypeCohort::sum_MACs() const {
  size_t sum = 0;
  for(size_t i = 0; i < get_n_sites(); i++) {
    sum += get_MAC(i);
  }
  return sum;
}

void haplotypeCohort::set_column(const vector<alleleValue>& values) {
  set_column(values, get_n_sites() - 1);
}

void haplotypeCohort::set_column(const vector<alleleValue>& values, size_t site) {
  haplotype_alleles_by_site_index[site] = values;
}

const siteIndex* haplotypeCohort::get_reference() const {
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
                                         reference->absolute_length(),
                                         reference->start_position() +
                                                  reference->absolute_length());
  haplotypeCohort* to_return = 
            new haplotypeCohort(number_of_haplotypes, new_reference);
  for(size_t i = 0; i < remaining_sites.size(); i++) {
    to_return->add_record(i);
    to_return->set_column(haplotype_alleles_by_site_index[remaining_sites[i]]);
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

vector<size_t> siteIndex::rand_positions(size_t N) const {
  return haploRandom::n_unique_uints(N, number_of_sites(), &site_index_to_position);
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
  vector<alleleValue> to_return(get_n_sites());
  size_t h_idx = rand_haplo_idx();
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  for(size_t i = 0; i < get_n_sites() - 1; i++) {
    binomial_distribution<size_t> recombiner(get_reference()->span_length_after(i) + 1, exp(log_recomb_probability));
    size_t n_recombinations = recombiner(generator);
    to_return[i] = haploRandom::mutate(allele_at(i, h_idx), log_mutation_probability);
    for(size_t i = 0; i < n_recombinations; i++) {
      h_idx = rand_haplo_idx(h_idx);
    }
  }
  to_return.back() = haploRandom::mutate(allele_at(get_n_sites() - 1, h_idx), log_mutation_probability);
  return to_return;
}

vector<alleleValue> haplotypeCohort::rand_desc_haplo(size_t generations, double log_recomb_probability, double log_mutation_probability) const {
  vector<vector<alleleValue> > old_gen;
  vector<vector<alleleValue> > new_gen;
  
  size_t gen_size = (size_t)pow(2, generations);
  new_gen.clear();
  for(size_t j = 0; j < gen_size; j++) {
    new_gen.push_back(rand_LS_haplo(log_recomb_probability, log_mutation_probability));
  }
  
  for(size_t i = 1; i < generations; i++) {
    gen_size = (size_t)pow(2, generations - i);
    old_gen = new_gen;
    new_gen.clear();
    for(size_t j = 0; j < gen_size; j++) {
      new_gen.push_back(get_reference()->make_child(old_gen[2*j], old_gen[2*j + 1], log_recomb_probability, log_mutation_probability));
    }
  }
  
  return(get_reference()->make_child(old_gen[0], old_gen[1], log_recomb_probability, log_mutation_probability));
}

vector<alleleValue> siteIndex::make_child(const vector<alleleValue>& parent_0, const vector<alleleValue>& parent_1, double log_recomb_probability, double log_mutation_probability) const {
  vector<alleleValue> to_return(number_of_sites());
  default_random_engine generator;
  generator.seed(chrono::system_clock::now().time_since_epoch().count());
  bernoulli_distribution hap_chooser(0.5);
  bool h_idx = hap_chooser(generator);
  for(size_t i = 0; i < number_of_sites() - 1; i++) {
    binomial_distribution<size_t> recombiner(span_length_after(i) + 1, exp(log_recomb_probability));
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
    if(draw >= (size_t)a) { draw++; }
    return (alleleValue)draw;
  }
  return a;
}