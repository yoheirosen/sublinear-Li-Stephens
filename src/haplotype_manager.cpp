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
              reference_sequence(referenceSequence(reference_bases)),
              cutoff_interval(thresholdInterval(penalties))
              {
  read_reference = string(read_bases);
  tree = new haplotypeStateTree(reference, penalties, cohort);
  end_position = start_position + read_reference.size() - 1;
  // TODO: handle case that read is not contained within reference
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

void haplotypeManager::set_cutoff_interval(double relative_threshold) {
  cutoff_interval = thresholdInterval(relative_threshold, penalties);
}

void haplotypeManager::build_entire_tree(double absolute_threshold) {
  initialize_tree();
  for(size_t i = 1; i < shared_sites(); i++) {
    build_next_level(absolute_threshold);
  }
  if(shared_sites() != 0) {
    extend_final_level(absolute_threshold);
  }
}


void haplotypeManager::build_entire_tree_interval(double cutoff) {
  initialize_tree();
  set_cutoff_interval(cutoff);
  cerr << "built root" << endl;
  cerr << "\tcurrent leaves " << current_leaves.size() << endl;
  for(size_t i = 1; i < shared_sites(); i++) {
    cerr << "building shared site " << i << " of "<< shared_sites() << "; ref index " << shared_index_to_ref_index(i) << endl;
    build_next_level_interval(0);
    cerr << "\tcurrent leaves " << current_leaves.size() << endl;
  }
  if(shared_sites() != 0) {
    extend_final_level(0);
  }
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
    // CASE 1: the first shared site is the first position in the read
    if(shared_sites() != 0) {
      if(get_shared_site_ref_position(0) == start_position) {
        // start_position is a shared site
        start_with_active_site(shared_index_to_ref_index(0));
        current_leaves = tree->root->get_unordered_children();
        return;
      }
    }
    // CASE 2: there exists sequence in the beginning of the read which is not a
    // shared site
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

void haplotypeManager::build_next_level_interval(double threshold) {
  if(last_level_built >= shared_sites() - 1) {
    // no more levels to build
    return;
  } else {
    size_t one_site_past_last_shared =    
              shared_index_to_ref_index(last_level_built) + 1;
    size_t current_site = shared_index_to_ref_index(last_level_built + 1);
    if(one_site_past_last_shared != current_site) {
      // extend all sites to current shared site
      fill_in_level_no_threshold(one_site_past_last_shared, current_site);
    }

    // copy past (smaller) vector to make space for new (larger) vector
    vector<haplotypeStateNode*> last_leaves = current_leaves;
    current_leaves.clear();

    vector<rowSet*> current_rows = get_rowSets_at_site(current_site);
    
    if(last_leaves.size() != 0) {
      branch_node(last_leaves[0], current_site, current_rows);
      cutoff_interval.set_new_site();
      cutoff_interval.check_for_new_bound(last_leaves[0]->get_unordered_children());
    }
    // thresholdInterval predictor(penalties);
    for(size_t i = 1; i < last_leaves.size(); i++) {
      haplotypeStateNode* n = last_leaves[i];
      branch_node_interval(n, current_site, current_rows);
      // branch_node_interval(n, current_site, current_rows, predictor);
    }
    for(size_t i = 0; i < last_leaves.size(); i++) {
      for(size_t j = 0; j < last_leaves[i]->number_of_children(); j++) {
        haplotypeStateNode* n = last_leaves[i]->get_child(j);
        if(!cutoff_interval.is_within_interval(n)) {
          n->mark_for_deletion();
        } else {
          current_leaves.push_back(n);
        }
      }
    }
    
    clear_rowSet_vector(current_rows);
    delete_marked_children(last_leaves);
    trim_back_abandoned_nodes(last_leaves);
  }
  last_level_built++;
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

    // copy past (smaller) vector to make space for new (larger) vector
    vector<haplotypeStateNode*> last_leaves = current_leaves;
    current_leaves.clear();

    vector<rowSet*> current_rows = get_rowSets_at_site(current_site);
    
    double likeliest_unrep_failure = threshold;
    for(size_t i = 0; i < last_leaves.size(); i++) {
      haplotypeStateNode* n = last_leaves[i];
      if(threshold == 0) {
        branch_node(n, current_site, current_rows);
        for(size_t j = 0; j < n->get_unordered_children().size(); j++) {
          current_leaves.push_back(n->get_unordered_children()[j]);
        }
      } else {
        if(!(n->is_marked_for_deletion())) {
          if(n->prefix_likelihood() >= threshold) {
            branch_node(n, current_site, current_rows, threshold, likeliest_unrep_failure);
            for(size_t j = 0; j < n->get_unordered_children().size(); j++) {
              haplotypeStateNode* n_child = n->get_unordered_children()[j];
              if(n_child->prefix_likelihood() > threshold) {
                current_leaves.push_back(n_child);
              } else {
                n_child->mark_for_deletion();
              }
            }
          } else {
            n->mark_for_deletion();
          }
        }
      }
    }
    clear_rowSet_vector(current_rows);
    delete_marked_children(last_leaves);
    trim_back_abandoned_nodes(last_leaves);
  }
  last_level_built++;
}

void haplotypeManager::delete_marked_children(vector<haplotypeStateNode*>& nodes) {
  for(size_t i = 0; i < nodes.size(); i++) {
    vector<haplotypeStateNode*> to_delete;
    for(size_t j = 0; j < nodes[i]->number_of_children(); j++) {
      if(nodes[i]->get_child(j)->is_marked_for_deletion()) {
        to_delete.push_back(nodes[i]->get_child(j));
      }
    }
    while(to_delete.size() > 0) {
      nodes[i]->remove_child(to_delete.back());
      to_delete.pop_back();
    }
  }
}

void haplotypeManager::trim_back_marked_nodes(vector<haplotypeStateNode*>& nodes) {
  vector<haplotypeStateNode*> to_delete;
  for(size_t i = 0; i < nodes.size(); i++) {
    if(nodes[i]->is_marked_for_deletion()) {
      to_delete.push_back(nodes[i]);
    }
  }
  nodes.clear();
  while(to_delete.size() > 0) {
    if(to_delete.back() != tree->root) {
      haplotypeStateNode* parent_of_deleted = to_delete.back()->get_parent();
      tree->remove_node(to_delete.back());
      to_delete.pop_back();
      if(parent_of_deleted->number_of_children() == 0) {
        to_delete.push_back(parent_of_deleted);
      }
    } else {
      to_delete.pop_back();
    }
  }
}

void haplotypeManager::trim_back_abandoned_nodes(vector<haplotypeStateNode*>& nodes) {
  for(size_t i = 0; i < nodes.size(); i++) {
    if(nodes[i]->number_of_children() == 0) {
      nodes[i]->mark_for_deletion();
    }
  }
  trim_back_marked_nodes(nodes);
}

void haplotypeManager::branch_node(haplotypeStateNode* n, size_t i) {
  fill_in_span_before(n, i);
  alleleValue a;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = n->add_child_copying_state(a);
    new_branch->state->extend_probability_at_site(i, a);
  }
  n->clear_state();
}

void haplotypeManager::branch_node(haplotypeStateNode* n, 
            size_t i, const vector<rowSet*>& rows, double threshold) {
  fill_in_span_before(n, i);
  alleleValue a;
  // if *anything* fails to pass threshold after extension, then this must also
  // happen for any allele not represented in the reference
  bool unrepresented_will_hit_threshold = false;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    if(threshold != 0) {
      if(!(unrepresented_will_hit_threshold && (cohort->number_matching(i, a) == 0))) {
        if(!will_hit_threshold(n, threshold, i, a)) {
          haplotypeStateNode* new_branch = n->add_child_copying_state(a);
          new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
          if(new_branch->prefix_likelihood() < threshold) {
            unrepresented_will_hit_threshold = true;
          }
        }
      }      
    } else {
      haplotypeStateNode* new_branch = n->add_child_copying_state(a);
      new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
    }
  }
  n->clear_state();
}

void haplotypeManager::branch_node_no_threshold(haplotypeStateNode* n, 
            size_t i, const vector<rowSet*>& rows) {
  fill_in_span_before(n, i);
  alleleValue a;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = n->add_child_copying_state(a);
    new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
  }
  n->clear_state();
}

void haplotypeManager::branch_node(haplotypeStateNode* n, 
            size_t i, const vector<rowSet*>& rows, double threshold, double& likeliest_unrep_failure) {
  fill_in_span_before(n, i);
  alleleValue a;
  // if *anything* fails to pass threshold after extension, then this must also
  // happen for any allele not represented in the reference
  bool unrepresented_will_hit_threshold = (n->prefix_likelihood() < likeliest_unrep_failure);
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    if(threshold != 0) {
      if(!(unrepresented_will_hit_threshold && (cohort->number_matching(i, a) == 0))) {
        if(!will_hit_threshold(n, threshold, i, a)) {
          haplotypeStateNode* new_branch = n->add_child_copying_state(a);
          new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
          if(new_branch->prefix_likelihood() < threshold) {
            likeliest_unrep_failure = new_branch->prefix_likelihood();
            unrepresented_will_hit_threshold = true;
          }
        }
      }      
    } else {
      haplotypeStateNode* new_branch = n->add_child_copying_state(a);
      new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
    }
  }
  n->clear_state();
}

void haplotypeManager::branch_node_interval(haplotypeStateNode* n, 
            size_t i, const vector<rowSet*>& rows) {
  fill_in_span_before(n, i);
  alleleValue a;

  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = n->add_child_copying_state(a);
    new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
    if(cutoff_interval.is_within_interval(new_branch)) {
      cutoff_interval.check_for_new_bound(new_branch);
    } else {
      new_branch->mark_for_deletion();
    }
  }
  n->clear_state();
}

void haplotypeManager::branch_node_interval(haplotypeStateNode* n, 
            size_t i, const vector<rowSet*>& rows, thresholdInterval& predictor) {
  fill_in_span_before(n, i);
  alleleValue a;
  // if *anything* fails to pass threshold after extension, then this must also
  // happen for any allele not represented in the reference
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = n->add_child_copying_state(a);
    if(!predictor.using_interval_cutoff() || predictor.is_within_interval(n)) {
      new_branch->state->extend_probability_at_site(*(rows[j]), get_cohort()->match_is_rare(i, a), a);
      if(cutoff_interval.is_within_interval(new_branch)) {
        cutoff_interval.check_for_new_bound(new_branch);
      } else {
        if(!predictor.using_interval_cutoff()) {
          predictor = thresholdInterval(n->prefix_likelihood(), penalties);
        } else {
          predictor.check_for_new_bound(n);
        }
        new_branch->mark_for_deletion();
      }
    }
  }
  n->clear_state();
}

void haplotypeManager::fill_in_level(double threshold,
        size_t start_site, size_t upper_bound_site) {
  if(threshold == 0) {
    fill_in_level_no_threshold(start_site, upper_bound_site);
  } else {
    fill_in_level_threshold(threshold, start_site, upper_bound_site);
  }
}

void haplotypeManager::fill_in_level_threshold(double threshold,
        size_t start_site, size_t upper_bound_site) {

  size_t p;
  alleleValue consensus_read_allele;
  haplotypeStateNode* n;
  
  for(size_t j = start_site; j < upper_bound_site; j++) {
    p = read_position(j);
    consensus_read_allele = char_to_allele(read_reference[p]);
    rowSet current_rowSet = get_cohort()->get_active_rowSet(j, consensus_read_allele);
    for(size_t i = 0; i < current_leaves.size(); i++) {
      if(current_leaves[i]->prefix_likelihood() < threshold) {
        current_leaves[i]->mark_for_deletion();
      } else if(!current_leaves[i]->is_marked_for_deletion()) {
        n = current_leaves[i];
        extend_node_at_site(n, j, consensus_read_allele, current_rowSet);
      }
    }
  }
}

void haplotypeManager::fill_in_level_no_threshold(size_t start_site, size_t upper_bound_site) {
  size_t p;
  alleleValue consensus_read_allele;
  haplotypeStateNode* n;
  
  for(size_t j = start_site; j < upper_bound_site; j++) {
    p = read_position(j);
    consensus_read_allele = char_to_allele(read_reference[p]);
    rowSet current_rowSet = get_cohort()->get_active_rowSet(j, consensus_read_allele);
    for(size_t i = 0; i < current_leaves.size(); i++) {
      n = current_leaves[i];
      extend_node_at_site(n, j, consensus_read_allele, current_rowSet);
    }
  }
}

void haplotypeManager::extend_final_level(double threshold) {
  if(threshold == 0) {
    extend_final_level_no_threshold();
  } else {
    extend_final_level_threshold(threshold);
  }
}

void haplotypeManager::extend_final_level_no_threshold() {
  if(ref_sites_after_shared_sites.back().size() != 0) {
    size_t past_last_shared = shared_index_to_ref_index(shared_sites() - 1) + 1;
    size_t past_last_ref = final_ref_site() + 1;
    fill_in_level(0, past_last_shared, past_last_ref);
  }
  if(final_span_after_last_ref_site() > 0) {
    haplotypeStateNode* n;
    for(size_t i = 0; i < current_leaves.size(); i++) {
      current_leaves[i]->state->extend_probability_at_span_after_anonymous(final_span_after_last_ref_site(), 0);
    }
  }
}

void haplotypeManager::extend_final_level_threshold(double threshold) {
  if(ref_sites_after_shared_sites.back().size() != 0) {
    size_t past_last_shared = shared_index_to_ref_index(shared_sites() - 1) + 1;
    size_t past_last_ref = final_ref_site() + 1;
    fill_in_level(threshold, past_last_shared, past_last_ref);
  }  

  if(final_span_after_last_ref_site() > 0) {
    vector<haplotypeStateNode*> last_leaves = current_leaves;
    current_leaves.clear();
    haplotypeStateNode* n;
    for(size_t i = 0; i < last_leaves.size(); i++) {
      n = last_leaves[i];
      if((!n->is_marked_for_deletion())) {
        n->state->extend_probability_at_span_after_anonymous(final_span_after_last_ref_site(), 0);
        if(n->prefix_likelihood() < threshold) {
          n->mark_for_deletion();
        } else {
          current_leaves.push_back(n);
        }
      }
    }
  }
}

void haplotypeManager::start_with_active_site(size_t i) {
  alleleValue a;
  for(size_t j = 0; j < 5; j++) {
    a = (alleleValue)j;
    haplotypeStateNode* new_branch = tree->root->add_child(a);
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

void haplotypeManager::extend_node_at_site(haplotypeStateNode* n, 
        size_t i, alleleValue a, const rowSet& row_set) {
  fill_in_span_before(n, i);
  n->state->extend_probability_at_site(row_set, cohort->match_is_rare(i, a), a);
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
  vector<size_t> unpruned_per_level;
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
      } else {
        level_prefix.push_back("");
      }
    }
    for(size_t i = 0; i < this_level.size(); i++) {
      if(!(this_level[i]->is_marked_for_deletion())) {   
        state_ID = tree->state_to_alleles(this_level[i]);
        for(size_t j = 0; j < state_ID.size(); j++) {
          cout << level_prefix[j] << "(" << allele_to_char(state_ID[j]) << ")";
        }
        // if(this_level[i]->is_marked_for_deletion()) {
          // cout << " : pruned" << endl;
        // } else {
          cout << " : " << this_level[i]->prefix_likelihood() << endl;
        // }
        if(this_level[i]->number_of_children() != 0) {
          temp_for_children = this_level[i]->get_unordered_children();
          for(size_t j = 0; j < temp_for_children.size(); j++) {
            next_level.push_back(temp_for_children[j]);
          }
        }
      }
    }
    level_depth++;
  }
  cout << total_nodes << " total nodes" << endl;
}

void haplotypeManager::print_tree_transitions() {
  vector<haplotypeStateNode*> next_level;
  vector<haplotypeStateNode*> this_level;
  vector<haplotypeStateNode*> temp_for_children;
  vector<alleleValue> state_ID;
  vector<string> level_prefix = {""};
  size_t level_depth = 0;
  size_t total_nodes = 0;
  vector<size_t> unpruned_per_level;
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
      } else {
        level_prefix.push_back("");
      }
    }
    for(size_t i = 0; i < this_level.size(); i++) {
      if(!(this_level[i]->is_marked_for_deletion())) {   
        state_ID = tree->state_to_alleles(this_level[i]);
        for(size_t j = 0; j < state_ID.size(); j++) {
          cout << level_prefix[j] << "(" << allele_to_char(state_ID[j]) << ")";
        }
        // if(this_level[i]->is_marked_for_deletion()) {
          // cout << " : pruned" << endl;
        // } else {
          cout << " : ";
          haplotypeStateNode* tracer = tree->root;
          for(size_t k = 0; k < state_ID.size(); k++) {
            tracer = tracer->get_child(state_ID[k]);
            double transition = tracer->prefix_likelihood() - tracer->max_prefix_likelihood(penalties);
            cout << transition << " ";
          }
          cout << endl;
        // }
        if(this_level[i]->number_of_children() != 0) {
          temp_for_children = this_level[i]->get_unordered_children();
          for(size_t j = 0; j < temp_for_children.size(); j++) {
            next_level.push_back(temp_for_children[j]);
          }
        }
      }
    }
    level_depth++;
  }
  cout << total_nodes << " total nodes" << endl;
}

vector<rowSet*> haplotypeManager::get_rowSets_at_site(size_t current_site) const {
	vector<rowSet*> to_return;
  for(size_t i = 0; i < 5; i++) {
    rowSet* to_add = new rowSet;
    *to_add = get_cohort()->get_active_rowSet(current_site, (alleleValue)i);
    to_return.push_back(to_add);
  }
  return to_return;
}

void haplotypeManager::clear_rowSet_vector(vector<rowSet*> row_sets) {
  for(size_t i = 0; i < 5; i++) {
    delete row_sets[i];
  }
}

bool haplotypeManager::will_hit_threshold(haplotypeStateNode* n, 
          double threshold, size_t site_index, alleleValue a) const {
  return ((n->prefix_likelihood() - threshold) < penalties->mu) &&
            (cohort->number_matching(site_index, a) == 0);
}

bool thresholdInterval::using_interval_cutoff() const {
  return using_interval;
}

thresholdInterval::thresholdInterval(const penaltySet* penalties) : penalties(penalties) {
  threshold = 0;
  using_interval = false;
}


thresholdInterval::thresholdInterval(double threshold,
          const penaltySet* penalties) : 
          threshold(threshold), penalties(penalties) {
  upper_bound = last_upper_bound - penalties->mu;
  if(threshold < 0) {
    using_interval = true;
  }
}

void thresholdInterval::set_new_site() {
  last_upper_bound = upper_bound;
  upper_bound = last_upper_bound + penalties->mu;
}

void thresholdInterval::check_for_new_bound(double test_bound) {
  if(test_bound > upper_bound) {
    upper_bound = test_bound;
  }
}

void thresholdInterval::check_for_new_bound(const haplotypeStateNode* test_bound) {
  if(test_bound->prefix_likelihood() > upper_bound) {
    upper_bound = test_bound->prefix_likelihood();
  }
}

void thresholdInterval::check_for_new_bound(const vector<double>& test_bounds) {
  upper_bound = test_bounds[0];
  for(size_t i = 1; i < test_bounds.size(); i++) {
    if(test_bounds[i] > upper_bound) {
      upper_bound = test_bounds[i];
    }
  }
}

void thresholdInterval::check_for_new_bound(const vector<haplotypeStateNode*>& test_bounds) {
  upper_bound = test_bounds[0]->prefix_likelihood();
  for(size_t i = 1; i < test_bounds.size(); i++) {
    if(test_bounds[i]->prefix_likelihood() > upper_bound) {
      upper_bound = test_bounds[i]->prefix_likelihood();
    }
  }
}

bool thresholdInterval::is_within_interval(double test_value) const {
  return test_value >= get_lower_bound();
}

bool thresholdInterval::is_within_interval(const haplotypeStateNode* test_value) const {
  return test_value->prefix_likelihood() >= get_lower_bound();
}

double thresholdInterval::get_upper_bound() const {
  return upper_bound;
}

double thresholdInterval::get_lower_bound() const {
  return upper_bound + threshold;
}

bool haplotypeManager::read_index_is_shared(size_t i) const {
  return read_site_is_shared[i];
}

bool haplotypeManager::read_matches(size_t i, alleleValue a) const {
  return reference_sequence.matches(i, a);
}

bool haplotypeManager::read_matches(size_t i, char a) const {
  return reference_sequence.matches(i, a);
}
