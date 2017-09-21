#include "haplotype_manager.hpp"
#include "set_of_extensions.hpp"

using namespace std;

haplotypeManager::haplotypeManager(
        const linearReferenceStructure* reference, const haplotypeCohort* cohort, 
              const penaltySet* penalties, const char* reference_bases,
        vector<size_t> site_positions_within_read,
        const char* read_bases, size_t start_reference_position) : 
        
        reference(reference), cohort(cohort), penalties(penalties),
        read_site_read_positions(site_positions_within_read),
              start_position(start_reference_position),
              reference_sequence(referenceSequence(reference_bases))
              {
  
  read_reference = string(read_bases);
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
}

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

size_t haplotypeManager::read_index_to_shared_index(size_t i) const {
  if(read_site_is_shared[i]) {
    return subsequence_indices[i];
  } else {
    return SIZE_MAX;
  }
}

size_t haplotypeManager::read_index_to_read_only_index(size_t i) const {
  if(!read_site_is_shared[i]) {
    return subsequence_indices[i];
  } else {
    return SIZE_MAX;
  }
}

size_t haplotypeManager::shared_index_to_read_index(size_t j) const {
  return shared_site_read_indices[j];
}

size_t haplotypeManager::shared_index_to_ref_index(size_t j) const {
  return get_ref_site_below_read_site(shared_index_to_read_index(j));
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

size_t haplotypeManager::final_ref_site_read_position() const {
  return read_position(get_ref_site_ref_position(final_ref_site()));
}

size_t haplotypeManager::final_shared_site_ref_index() const {
  return shared_index_to_ref_index(shared_sites() - 1);
}

size_t haplotypeManager::final_ref_site() const {
  if(!ref_sites) {
    return SIZE_MAX;
  } else {
    if(shared_sites()) {
      if(ref_sites_after_shared_sites.back().size() == 0) {
        return final_shared_site_ref_index();
      } else {
        return ref_sites_after_shared_sites.back().back().site_index;
      }
    } else {
      return ref_sites_in_initial_span.back().site_index;
    }
  }
}

size_t haplotypeManager::final_read_site_read_index() const {
  return read_sites() - 1;
}

size_t haplotypeManager::final_read_site_read_position() const {
  return get_read_site_read_position(read_sites() - 1);
}

size_t haplotypeManager::final_shared_site_read_index() const {
  return shared_index_to_read_index(shared_sites() - 1);
}

size_t haplotypeManager::final_shared_site_read_position() const {
  return get_read_site_read_position(final_shared_site_read_index());
}

size_t haplotypeManager::final_span_after_last_ref_site() const {
  return end_position - get_ref_site_ref_position(final_ref_site());
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

void haplotypeManager::find_ref_only_sites_and_alleles() {
  if(!contains_ref_sites()) {
    return;
  } else {
    size_t lower_ref_index;
    size_t upper_ref_index;
    
    // initial span
    lower_ref_index = reference->find_site_above(start_position);
    upper_ref_index = shared_index_to_ref_index(0);
    vector<alleleAtSite> to_add;
    for(size_t i = lower_ref_index; i < upper_ref_index; i++) {
      to_add.push_back(alleleAtSite(i, 
              char_to_allele(read_reference.at(read_position(
                      reference->get_position(i))))));
    }
    ref_sites_in_initial_span = to_add;
    
    // spans following read sites i to 1-before-end
    for(size_t i = 0; i < shared_sites() - 1; i++) {
      lower_ref_index = shared_index_to_ref_index(i) + 1;
      upper_ref_index = shared_index_to_ref_index(i + 1);
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
      lower_ref_index = shared_index_to_ref_index(shared_sites() - 1);
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
      if(!reference_sequence.matches(p, 
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
        if(!reference_sequence.matches(p, 
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
        if(!reference_sequence.matches(p, 
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
      if(!reference_sequence.matches(p, 
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
          if(!reference_sequence.matches(p, 
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
        if(!reference_sequence.matches(p, 
                char_to_allele(read_reference.at(read_position(p))))) {
          running_count++;
        }
      }
      invariant_penalties_by_ref_site.push_back(running_count);
    }
  }
}

void haplotypeManager::check_for_ref_sites() {
  ref_sites = (reference->find_site_above(start_position) != 
          reference->find_site_above(end_position));
}

size_t haplotypeManager::get_shared_site_ref_position(size_t j) const {
  return get_ref_site_ref_position(shared_index_to_ref_index(0));
}

size_t haplotypeManager::get_shared_site_read_position(size_t j) const {
  return get_read_site_read_position(shared_index_to_read_index(j));
}

size_t haplotypeManager::get_ref_site_ref_position(size_t j) const {
  return reference->get_position(j);
}

void haplotypeManager::initialize_tree() {
  last_level_built = 0;
  if(!contains_ref_sites()) {
    // there are zero reference sites in all of the read. Clearly no shared 
    // sites either

    // initialize tree which consists in its entirety of a single root node
    // whose state is given by a site-less span of length length()
    size_t read_length = length();
    start_with_span(read_length);
    current_leaves = {tree->root};
  } else {
    if(shared_sites() != 0) {
      if(get_shared_site_ref_position(0) == start_position) {
        // start_position is a shared site
        start_with_active_site(shared_index_to_ref_index(0));
        current_leaves = tree->root->get_unordered_children();
        return;
      }
    }
    if(ref_sites_in_initial_span.size() == 0) {
      // it's just a regular span from the point of view of the reference
      size_t initial_span_length = 
                get_shared_site_ref_position(0) == start_position;
      start_with_span(initial_span_length);
      current_leaves = {tree->root};
    } else {
      size_t first_ref = ref_sites_in_initial_span[0].site_index;
      size_t first_ref_pos = get_ref_site_ref_position(first_ref);
      alleleValue first_ref_allele = ref_sites_in_initial_span[0].allele;
      size_t next_ref;
      alleleValue next_ref_allele;
      
      // handle first reference site among those preceding the first shared site
      if(first_ref_pos == start_position) {
        start_with_inactive_site(first_ref, first_ref_allele);
      } else {
        size_t initial_span_length = first_ref_pos - start_position;
        // cout << initial_span_length << " initial span" << endl; // OK
        start_with_span(initial_span_length);
        extend_node_at_site(tree->root, first_ref, first_ref_allele);
      }
      
      // handle subsequent reference sites
      for(size_t i = 1; i < ref_sites_in_initial_span.size(); i++) {
        next_ref = ref_sites_in_initial_span[i].site_index;
        next_ref_allele = ref_sites_in_initial_span[i].allele;
        extend_node_at_site(tree->root, next_ref, next_ref_allele);
      }
      
      if(shared_sites() != 0) {
        branch_node(tree->root, shared_index_to_ref_index(0));
        current_leaves.clear();
        current_leaves = tree->root->get_unordered_children();
      } else {
        current_leaves = {tree->root};
      }
    }
  }
}

void haplotypeManager::build_next_level(double threshold) {
  if(last_level_built >= shared_sites() - 1) {
    // no more levels to build
    return;
  } else {
    size_t one_site_past_last_shared =    
              shared_index_to_ref_index(last_level_built) + 1;
    size_t current_site = shared_index_to_ref_index(last_level_built + 1);
    if(one_site_past_last_shared != current_site) {
      // extend all sites to current shared site
      fill_in_level(threshold, one_site_past_last_shared, current_site);
    }
    // control scope to avoid double-delete
    {
      // copy past (smaller) vector to make space for new (larger) vector
      vector<haplotypeStateNode*> last_leaves = current_leaves;
      current_leaves.clear();
      for(size_t i = 0; i < last_leaves.size(); i++) {
        haplotypeStateNode* n = last_leaves[i];
        // we may have deleted leaves from the previous level if they were found
        // to score below the threshold during the extension over the intervening
        // "spans." Need to check for this
        branch_node(n, current_site);
        for(size_t j = 0; j < n->get_unordered_children().size(); j++) {
          current_leaves.push_back(n->get_unordered_children()[j]);
        }
      }
    }
    last_level_built++;
  }
}

// TODO pre-estimate killed nodes
// make sure you account for read-ref mismatches
void haplotypeManager::fill_in_level(double threshold,
        size_t start_site, size_t upper_bound_site) {

  size_t one_site_past_last_shared =    
            shared_index_to_ref_index(last_level_built) + 1;
  size_t current_site = shared_index_to_ref_index(last_level_built + 1);
  size_t p;
  alleleValue consensus_read_allele;
  haplotypeStateNode* n;
  
  for(size_t j = start_site; j < upper_bound_site; j++) {
    p = read_position(j);
    consensus_read_allele = char_to_allele(read_reference[p]);
    if(threshold != 0) {
      //TODO
      for(size_t i = 0; i < current_leaves.size(); i++) {
        if(current_leaves[i]->prefix_likelihood() < 0) {
          
        }
      }
    } else {
      for(size_t i = 0; i < current_leaves.size(); i++) {
        n = current_leaves[i];
        extend_node_at_site(n, j, consensus_read_allele);
      }
    }
  }
}

void haplotypeManager::extend_final_level(double threshold) {
  if(ref_sites_after_shared_sites.back().size() != 0) {
    size_t past_last_shared = shared_index_to_ref_index(shared_sites() - 1) + 1;
    size_t past_last_ref = final_ref_site() + 1;
    fill_in_level(threshold, past_last_shared, past_last_ref);
  }
  vector<haplotypeStateNode*> temp = current_leaves;
  current_leaves.clear();
  haplotypeStateNode* n;
  if(final_span_after_last_ref_site() > 0) {
    for(size_t i = 0; i < temp.size(); i++) {
      n = temp[i];
      if(n != nullptr) {
        n->state->extend_probability_at_span_after_anonymous(final_span_after_last_ref_site(), 0);
        if(threshold != 0) {
          if(n->state->prefix_likelihood() < threshold) {
            tree->remove_node_and_unshared_ancestors(n);
          } else {
            current_leaves.push_back(n);
          }
        } else {
          current_leaves.push_back(n);
        }
      }
    }
  }
}

void haplotypeManager::build_entire_tree(double threshold) {
  initialize_tree();
  for(size_t i = 1; i < shared_sites(); i++) {
    build_next_level(threshold);
  }
  if(shared_sites() != 0) {
    extend_final_level(threshold);
  }
}

void haplotypeManager::start_with_active_site(size_t i) {
  alleleValue a;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = tree->root->add_child(a);
    // start_sites.push_back(new_branch);
    new_branch->state = new haplotypeMatrix(reference, penalties, cohort);
    new_branch->state->initialize_probability_at_site(i, a);
  }
}

void haplotypeManager::start_with_inactive_site(size_t i, alleleValue a) {
  tree->root->state = new haplotypeMatrix(reference, penalties, cohort);
  tree->root->state->initialize_probability_at_site(i, a);
}

void haplotypeManager::start_with_span(size_t length) {
  tree->root->state = new haplotypeMatrix(reference, penalties, cohort);
  tree->root->state->initialize_probability_at_span(length, 0);
}

void haplotypeManager::fill_in_span_before(haplotypeStateNode* n, size_t i) {
  if(!(n->state->last_extended_is_span()) && reference->has_span_before(i)) {
    n->state->extend_probability_at_span_after(i-1, 0);
  }
}

void haplotypeManager::extend_node_at_site(haplotypeStateNode* n, 
        size_t i, alleleValue a) {
  fill_in_span_before(n, i);
  n->state->extend_probability_at_site(i, a);
}

void haplotypeManager::branch_node(haplotypeStateNode* n, 
        size_t i) {
  fill_in_span_before(n, i);
  alleleValue a;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = n->add_child_copying_state(a);
    new_branch->state->extend_probability_at_site(i, a);
  }
  n->clear_state();
}

vector<haplotypeStateNode*> haplotypeManager::get_current_leaves() const {
  return current_leaves;
}

void haplotypeManager::print_tree() {
  vector<haplotypeStateNode*> next_level;
  vector<haplotypeStateNode*> this_level;
  vector<haplotypeStateNode*> temp_for_children;
  vector<alleleValue> state_ID;
  vector<string> level_prefix = {""};
  size_t level_depth = 0;
  size_t total_nodes = 0;
  cout << "root : " << tree->root->prefix_likelihood() << endl;
  if(tree->root->number_of_children() != 0) {
    next_level = tree->root->get_unordered_children();
  }
  while(next_level.size() != 0) {
    this_level = next_level;
    next_level.clear();
    total_nodes += this_level.size();
    if(level_depth == 0) {
      if(get_shared_site_read_position(0) != 0) {
        level_prefix[0] = read_reference.substr(0, 
                    get_shared_site_read_position(0));
      }
    } else {
      if(get_shared_site_read_position(level_depth) - 
                  get_shared_site_read_position(level_depth - 1) > 1) {
        size_t bdd_1 = get_shared_site_read_position(level_depth - 1) + 1;
        size_t len = get_shared_site_read_position(level_depth) - bdd_1;
        level_prefix.push_back(read_reference.substr(bdd_1, len));
      }
    }
    for(size_t i = 0; i < this_level.size(); i++) {   
      state_ID = tree->state_to_alleles(this_level[i]);
      for(size_t j = 0; j < state_ID.size(); j++) {
        cout << level_prefix[j] << "(" << allele_to_char(state_ID[j]) << ")";
      }
      cout << " : " << this_level[i]->prefix_likelihood() << endl;
      if(this_level[i]->number_of_children() != 0) {
        temp_for_children = this_level[i]->get_unordered_children();
        for(size_t j = 0; j < temp_for_children.size(); j++) {
          next_level.push_back(temp_for_children[j]);
        }
      }
    }
    level_depth++;
  }
  cout << total_nodes << " total nodes" << endl;
}