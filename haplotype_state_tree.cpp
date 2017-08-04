#include "haplotype_state_tree.hpp"

using namespace std;

haplotypeStateTree::haplotypeStateTree() {
  root = new haplotypeStateNode();
}

haplotypeStateTree::haplotypeStateTree(linearReferenceStructure* reference, 
            penaltySet* pen, haplotypeCohort* cohort) :
            ref(reference), pen(pen), cohort(cohort) {
  root = new haplotypeStateNode();
}

haplotypeStateTree::~haplotypeStateTree() {
  delete root;
}

haplotypeStateNode* haplotypeStateTree::alleles_to_state(
            vector<alleleValue> identifiers) {
  haplotypeStateNode* current_node = root;
  for(size_t i = 0; i < identifiers.size(); i++) {
    if(current_node->is_leaf()) {
      current_node == nullptr;
      break;
    }
    current_node = current_node->get_child(identifiers[i]);
    if(current_node == nullptr) {
      break;
    }
  }
  return current_node;
}

vector<alleleValue> haplotypeStateTree::state_to_alleles(
            haplotypeStateNode* state_node) {
  haplotypeStateNode* current_node = state_node;
  vector<alleleValue> backwards_allele_vector;
  while(current_node != root) {
    backwards_allele_vector.push_back(current_node->identifying_allele());
    current_node = current_node->get_parent();
  }
  backwards_allele_vector.push_back(root->identifying_allele());
  vector<alleleValue> to_return;
  size_t upper = backwards_allele_vector.size() - 1;
  for(size_t i = 0; i < backwards_allele_vector.size(); i++) {
    to_return.push_back(backwards_allele_vector[upper - i]);
  }
  return to_return;
}

void haplotypeStateTree::remove_node(haplotypeStateNode* n) {
  if(n == root) {
    delete root;
  } else {
    haplotypeStateNode* parent = n->get_parent();
    // handles both deletion from parent's child-vector as well as 
    parent->remove_child(n);
  }
}

void haplotypeStateTree::remove_node_and_unshared_ancestors(
            haplotypeStateNode* n) {
  if(!n->is_leaf()) {
    remove_node(n);
  } else {
    haplotypeStateNode* current_node = n;
    haplotypeStateNode* next_node;
    while(current_node->is_leaf()) {
      next_node = current_node->get_parent();
      // remove from children of parent
      delete current_node;
      current_node = next_node;
    }
    current_node->remove_child_from_childvector();
  }
}

void haplotypeStateTree::remove_unlikely_children(haplotypeStateNode* n, 
            double threshold) {
  for(size_t i = 0; i < n->number_of_children(); i++) {
    if(n->get_child(i)->prefix_likelihood() < threshold) {
      n->remove_child(n->get_child(i));
    }
  }
}

size_t haplotypeStateTree::length_to_node(haplotypeStateNode* n) const {
  return ref->get_position(n->get_end_site_index()) - initial_position;
}

void haplotypeStateTree::start_with_active_site(size_t i) {
  segregating_sites.push_back(i);
  for(size_t j = 0; j < reference->get_alleles_at_site().size(); j++) {
    haplotypeStateNode* new_branch = 
            root->add_child(alleleAtSite(i, 
                    reference->get_alleles_at_site()[j]));
  }
}

void haplotypeStateTree::start_with_inactive_site(size_t i, alleleValue a) {
  root->state = new haplotypeMatrix(ref, pen, cohort);
  root->state->initialize_probability_at_site(i, a);
}

void haplotypeStateTree::start_with_span(size_t length) {
  root->state = new haplotypeMatrix(ref, pen, cohort);
  root->state->initialize_probability_at_span(length, 0);
}

void haplotypeStateTree::fill_in_span_before(haplotypeStateNode* n, size_t i) {
  if(!(n->state->last_extended_is_span()) && n->state->has_span_before(i)) {
    n->state->extend_probability_at_span_after(i-1, 0);
  }
}

void haplotypeStateTree::extend_node_by_allele_at_site(haplotypeStateNode* n, 
        size_t i, alleleValue a) {
  fill_in_span_before(n, i);
  n->state->extend_probability_at_site(i, a);
}

void haplotypeStateTree::branch_node_by_alleles_at_site(haplotypeStateNode* n, 
        size_t i, double threshold) {
  fill_in_span_before(n, i);
  for(size_t j = 0; j < reference->get_alleles_at_site().size(); j++) {
    haplotypeStateNode* new_branch = 
            n->add_child_copying_state(alleleAtSite(i, 
                    reference->get_alleles_at_site()[j]));
  }
  delete n->state;
  if(threshold != 0) {
    remove_unlikely_children(n, threshold);
  }
}