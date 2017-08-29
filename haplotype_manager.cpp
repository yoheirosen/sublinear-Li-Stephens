#include "haplotype_manager.hpp"

using namespace std;

size_t haplotypeManager::length() const {
  return read_reference.length();
}

size_t haplotypeManager::read_sites() const {
  return read_site_read_positions.size();
}

size_t haplotypeManager::shared_sites() const {
  return shared_site_read_indices.size();
}

size_t haplotypeManager::ref_position(size_t p) const {
  return p + start_position;
}

size_t haplotypeManager::read_position(size_t p) const {
  if(p >= start_position && p <= end_position) {
    return p - start_position;
  } else {
    return SIZE_MAX;
  }
}

size_t haplotypeManager::get_read_site_read_position(size_t i) const {
  return read_site_read_positions[i];
}

size_t haplotypeManager::get_read_site_ref_position(size_t i) const {
  return ref_position(read_site_read_positions[i]);
}

void haplotypeManager::find_shared_sites() {
  for(size_t i = 0; i < read_site_read_positions.size(); i++) {
    read_site_is_shared.push_back(
            reference->is_site(
                  ref_position(read_site_read_positions[i])));
    if(read_site_is_shared[i]) {
      shared_site_read_indices.push_back(i);
    }
  }
}

void haplotypeManager::find_ref_sites_below_read_sites() {
  for(size_t i = 0; i < read_site_read_positions.size(); i++) {
    ref_site_below_read_site.push_back(
            reference->find_site_below(
                    ref_position(read_site_read_positions[i])));
  }
}

void haplotypeManager::build_subsequence_indices() {
  size_t next_read_only = 0;
  size_t next_shared = 0;
  for(size_t i = 0; i < read_site_read_positions.size(); i++) {
    if(read_site_is_shared[i]) {
      subsequence_indices.push_back(next_shared);
      next_shared++;
    } else {
      subsequence_indices.push_back(next_read_only);
      next_read_only++;
    }
  }
}

size_t haplotypeManager::index_among_shared_sites(size_t i) const {
  if(read_site_is_shared[i]) {
    return subsequence_indices[i];
  } else {
    return SIZE_MAX;
  }
}

size_t haplotypeManager::index_among_read_only_sites(size_t i) const {
  if(!read_site_is_shared[i]) {
    return subsequence_indices[i];
  } else {
    return SIZE_MAX;
  }
}

size_t haplotypeManager::get_shared_site_read_index(size_t j) const {
  return shared_site_read_indices[j];
}

size_t haplotypeManager::get_shared_site_ref_index(size_t j) const {
  return get_ref_site_below_read_site(get_shared_site_read_index(j));
}

size_t haplotypeManager::get_ref_site_below_read_site(size_t i) const {
  return ref_site_below_read_site[i];
}

double haplotypeManager::invariant_penalty_at_read_site(size_t i) const {
  if(read_reference.size() == 0) {
    return 0;
  } else {
    return (penalties->mu)*invariant_penalties_by_read_site[i];
  }
}

double haplotypeManager::invariant_penalty_at_ref_site(size_t i) const {
  if(read_reference.size() == 0) {
    return 0;
  } else {
    return (penalties->mu)*invariant_penalties_by_ref_site[i];
  }
}

bool haplotypeManager::contains_shared_sites() const {
  return (shared_site_read_indices.size() != 0);
}

bool haplotypeManager::contains_ref_sites() const {
  return ref_sites;
}

bool haplotypeManager::contains_read_only_sites() const {
  return (shared_site_read_indices.size() != read_site_read_positions.size());
}

void haplotypeManager::check_for_ref_sites() {
  ref_sites = (reference->find_site_above(start_position) != 
          reference->find_site_above(end_position));
}

void haplotypeManager::find_ref_only_sites_and_alleles() {
  if(!contains_ref_sites()) {
    return;
  } else {
    size_t lower_ref_index;
    size_t upper_ref_index;
    
    // initial span
    lower_ref_index = reference->find_site_above(start_position);
    upper_ref_index = get_shared_site_ref_index(0);
    vector<alleleAtSite> to_add;
    for(size_t i = lower_ref_index; i < upper_ref_index; i++) {
      to_add.push_back(alleleAtSite(i, 
              char_to_allele(read_reference.at(read_position(
                      reference->get_position(i))))));
    }
    ref_sites_in_initial_span = to_add;
    
    // spans following read sites i to 1-before-end
    for(size_t i = 0; i < shared_sites() - 1; i++) {
      lower_ref_index = get_shared_site_ref_index(i) + 1;
      upper_ref_index = get_shared_site_ref_index(i + 1);
      to_add.clear();
      for(size_t j = lower_ref_index; j < upper_ref_index; j++) {
        to_add.push_back(alleleAtSite(j, 
                char_to_allele(read_reference.at(read_position(
                        reference->get_position(j))))));
      }
      ref_sites_after_shared_sites.push_back(to_add);
    }
    
    // terminal span
    if(shared_sites() != 0) {
      lower_ref_index = get_shared_site_ref_index(shared_sites() - 1);
      upper_ref_index = reference->find_site_above(end_position);
      to_add.clear();
      for(size_t j = lower_ref_index; j < upper_ref_index; j++) {
        to_add.push_back(alleleAtSite(j, 
                char_to_allele(read_reference.at(read_position(
                        reference->get_position(j))))));
      }
      ref_sites_after_shared_sites.push_back(to_add);
    }
  }
}

void haplotypeManager::count_invariant_penalties() {
  if(read_reference.size() == 0) {
    // read reference is the same as the reference structure reference--the only
    // possible deviations are at the read-sites
    return;
  } else {
    size_t running_count = 0;
    size_t count_from;
    size_t count_until;
    
    // count by read site
    // initial span
    count_from = start_position;    
    if(read_site_read_positions.size() != 0) {
      count_until = get_read_site_ref_position(0);
    } else {
      count_until = end_position + 1;
    }
    for(size_t p = count_from; p < count_until; p++) {
      if(!reference_sequence->matches(p, 
              char_to_allele(read_reference.at(read_position(p))))) {
        running_count++;
      }
    }
    invariant_penalties_by_read_site.push_back(running_count);
    
    // spans following read sites i to 1-before-end
    for(size_t i = 0; i < read_site_read_positions.size() - 1; i++) {
      count_from = get_read_site_ref_position(i);
      count_until = get_read_site_ref_position(i + 1);
      for(size_t p = count_from; p < count_until; p++) {
        if(!reference_sequence->matches(p, 
                char_to_allele(read_reference.at(read_position(p))))) {
          running_count++;
        }
      }
      invariant_penalties_by_read_site.push_back(running_count);
    }
    
    // terminal span
    if(read_site_read_positions.size() != 0) {
      count_from = get_read_site_ref_position(read_sites() - 1);
      count_until = end_position + 1;
      for(size_t p = count_from; p < count_until; p++) {
        if(!reference_sequence->matches(p, 
                char_to_allele(read_reference.at(read_position(p))))) {
          running_count++;
        }
      }
      invariant_penalties_by_read_site.push_back(running_count);
    }
    
    // count by ref site
    running_count = 0;
    
    count_from = start_position;    
    if(contains_ref_sites()) {
      count_until = reference->get_position(
              reference->find_site_above(start_position));
    } else {
      count_until = end_position + 1;
    }
    for(size_t p = count_from; p < count_until; p++) {
      if(!reference_sequence->matches(p, 
              char_to_allele(read_reference.at(read_position(p))))) {
        running_count++;
      }
    }
    invariant_penalties_by_ref_site.push_back(running_count);
    
    if(contains_ref_sites()) {
      // spans following ref sites i to 1-before-end
      for(size_t i = reference->find_site_above(start_position) + 1;
              i < reference->find_site_below(end_position) - 1; i++) {
        count_from = reference->get_position(i) + 1;
        count_until = reference->get_position(i + 1) - 1;
        for(size_t p = count_from; p < count_until; p++) {
          if(!reference_sequence->matches(p, 
                  char_to_allele(read_reference.at(read_position(p))))) {
            running_count++;
          }
        }
        invariant_penalties_by_ref_site.push_back(running_count);
      }
      
      count_from = 
             reference->get_position(
                    reference->find_site_below(end_position)) + 1;
      count_until = end_position + 1;
      for(size_t p = count_from; p < count_until; p++) {
        if(!reference_sequence->matches(p, 
                char_to_allele(read_reference.at(read_position(p))))) {
          running_count++;
        }
      }
      invariant_penalties_by_ref_site.push_back(running_count);
    }
  }
}

haplotypeManager::haplotypeManager(
        const linearReferenceStructure* reference, const haplotypeCohort* cohort, 
              const penaltySet* penalties, const char* reference_bases,
        vector<size_t> site_positions_within_read,
        const char* read_bases, size_t start_reference_position) : 
        
        reference(reference), cohort(cohort), penalties(penalties),
        read_site_read_positions(site_positions_within_read),
        read_reference(read_reference), 
              start_position(start_reference_position) {
  
  read_reference = string(read_bases);
  reference_sequence = new referenceSequence(reference_bases);
  tree = new haplotypeStateTree(reference, penalties, cohort);
  end_position = start_position + read_reference.size() - 1;
  find_ref_sites_below_read_sites();
  find_shared_sites();
  check_for_ref_sites();
  build_subsequence_indices();
  count_invariant_penalties();
  find_ref_only_sites_and_alleles();
}

haplotypeManager::~haplotypeManager() {
  delete tree;
  delete reference_sequence;
}

void haplotypeManager::initialize_tree() {
  if(!contains_ref_sites()) {
    // there are zero reference sites in all of the read. Clearly no shared 
    // sites either
    
    // initialize tree which consists in its entirety of a single root node
    // whose state is given by a site-less span of length length()
    
    size_t augmentations = invariant_penalties_by_read_site[0];
    size_t read_length = length();
    tree->start_with_span(read_length);
  } else {
    if(reference->get_position(get_shared_site_ref_index(0)) ==
            start_position) {
      return;
    } else {
      if(ref_sites_in_initial_span.size() == 0) {
        // it's just a regular span
        tree->start_with_span(
                reference->get_position(get_shared_site_ref_index(0)) -
                start_position);
      } else {
        // there are some ref sites in here
        if(reference->get_position(ref_sites_in_initial_span[0].site_index) ==
                start_position) {
          tree->start_with_inactive_site(
                  ref_sites_in_initial_span[0].site_index,
                  ref_sites_in_initial_span[0].allele);
        } else {
          size_t left_span_length =
                  reference->get_position(
                          ref_sites_in_initial_span[0].site_index) -
                  start_position;
          tree->start_with_span(left_span_length);
          tree->extend_node_by_allele_at_site(tree->root,
                  ref_sites_in_initial_span[0].site_index,
                  ref_sites_in_initial_span[0].allele);
        }
        for(size_t i = 1; i < ref_sites_in_initial_span.size(); i++) {
          tree->extend_node_by_allele_at_site(tree->root,
                  ref_sites_in_initial_span[i].site_index,
                  ref_sites_in_initial_span[i].allele);
        }
      }
    }
  }
}

void haplotypeManager::build_first_level(double threshold) {
  initialize_tree();
  if(shared_sites() > 0) {
    tree->branch_node_by_alleles_at_site(tree->root,
            get_shared_site_ref_index(0), threshold);
    nodes_at_last_level_built = tree->root->get_unordered_children();
    last_level_built = 0;
  }
}

void haplotypeManager::build_next_level(double threshold) {
  if(last_level_built >= shared_sites() - 1) {
    return;
  } else {
    // step all states forward
    // TODO: confirm this things meets neccessary gurantees
    fill_in_level(threshold, get_shared_site_ref_index(last_level_built) + 1,
            get_shared_site_ref_index(last_level_built + 1));
    // TODO reserve memory for this massive thing! And be smarter with copying
    vector<haplotypeStateNode*> temp;
    double current_threshold = threshold
            - invariant_penalty_at_read_site(
                      get_shared_site_read_index(
                                last_level_built + 1));
    for(size_t i = 0; i < nodes_at_last_level_built.size(); i++) {
      haplotypeStateNode* n = nodes_at_last_level_built[i];
      if(n != nullptr) {
        tree->branch_node_by_alleles_at_site(n, 
                  get_shared_site_ref_index(last_level_built + 1), 
                          current_threshold);
        if(n->number_of_children() == 0) {
          tree->remove_node_and_unshared_ancestors(n);
        } else {
          for(size_t j = 0; j < n->get_unordered_children().size(); j++) {
            temp.push_back(n->get_unordered_children()[j]);
          }
        }
      }
    }
    last_level_built++;
  }
}

void haplotypeManager::fill_in_level(double threshold,
        size_t start_site, size_t upper_bound_site) {
  double highest_failure = 0;
  for(size_t i = 0; i < nodes_at_last_level_built.size(); i++) {
    haplotypeStateNode* n = nodes_at_last_level_built[i];
    double starting_level = n->prefix_likelihood();
    if(starting_level < highest_failure - penalties->rho
            && highest_failure != 0) {
      tree->remove_node_and_unshared_ancestors(n);
      nodes_at_last_level_built[i] = nullptr;
    } else {
      // TODO: pre-computation bounds per-allele
      for(size_t j = start_site; j < upper_bound_site; j++) {
        alleleValue a =
                (ref_sites_after_shared_sites[last_level_built + 1][
                j - get_shared_site_ref_index(last_level_built) - 1]).allele;
        tree->extend_node_by_allele_at_site(n, j, a);
        if(n->prefix_likelihood() + invariant_penalty_at_ref_site(j)
                < threshold) {
          tree->remove_node_and_unshared_ancestors(n);
          nodes_at_last_level_built[i] = nullptr;
          if(highest_failure == 0) {
            highest_failure = starting_level;
          } else {
            if(starting_level > highest_failure) {
              highest_failure = starting_level;
            }
          }
          break;
        }
      }
    }
  }
}

void haplotypeManager::extend_final_level(double threshold) {
  fill_in_level(threshold,
          get_shared_site_ref_index(shared_sites() - 1),
          ref_sites_after_shared_sites.back().back().site_index + 1);
  vector<haplotypeStateNode*> temp = nodes_at_last_level_built;
  nodes_at_last_level_built.clear();
  size_t terminal_span = end_position;
  if(terminal_span != 0) {
    for(size_t i = 0; i < temp.size(); i++) {
      haplotypeStateNode* n = temp[i];
      if(n != nullptr) {
        // TODO: variable length final span function
      }
    }
  }
}

void haplotypeManager::build_entire_tree(double threshold) {
  build_first_level(threshold);
  for(size_t i = 1; i < shared_sites(); i++) {
    build_next_level(threshold);
  }
  if(shared_sites() != 0) {
    extend_final_level(threshold);
  }
}

size_t haplotypeManager::levels_built() const {
  return last_level_built;
}

bool haplotypeManager::all_levels_built() const {
  return last_level_built == length();
}

const haplotypeStateTree* haplotypeManager::get_tree() const {
  return tree;
}

const linearReferenceStructure* haplotypeManager::get_reference() const {
  return reference;
}

const haplotypeCohort* haplotypeManager::get_cohort() const {
  return cohort;
}

const penaltySet* haplotypeManager::get_penalties() const {
  return penalties;
}

haplotypeStateNode* haplotypeManager::find_node_by_prefix(string& prefix) {
  vector<alleleValue> prefix_alleles;
  for(size_t i = 0; i < prefix.length(); i++) {
    prefix_alleles.push_back(char_to_allele(prefix[i]));
  }
  haplotypeStateNode* result = tree->alleles_to_state(prefix_alleles);
}