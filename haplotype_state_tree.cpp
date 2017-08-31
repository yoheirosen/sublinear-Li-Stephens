#include "haplotype_state_tree.hpp"

using namespace std;

haplotypeStateTree::haplotypeStateTree() {
  root = new haplotypeStateNode();
}

haplotypeStateTree::haplotypeStateTree(const linearReferenceStructure* reference, 
            const penaltySet* penalties, const haplotypeCohort* cohort) :
            reference(reference), penalties(penalties), cohort(cohort) {
  root = new haplotypeStateNode();
}

haplotypeStateTree::~haplotypeStateTree() {
  delete root;
}

// TODO test []
void haplotypeStateTree::remove_node(haplotypeStateNode* n) {
  if(n == root) {
    delete root;
  } else {
    haplotypeStateNode* parent = n->get_parent();
    // handles both deletion from parent's child-vector as well as 
    parent->remove_child(n);
  }
}

// TODO test []
void haplotypeStateTree::remove_node_and_unshared_ancestors(
            haplotypeStateNode* n) {
  if(n == root) {
    remove_node(n);
  } else {
    haplotypeStateNode* current_node = n->get_parent();
    haplotypeStateNode* next_node;
    remove_node(n);
    while(current_node->is_leaf()) {
      if(current_node == root) {
        remove_node(n);
        break;
      } else {
        next_node = current_node->get_parent();
        // remove from children of parent
        next_node->remove_child(current_node);
        current_node = next_node;
      }
    }
  }
}

// TODO test []
void haplotypeStateTree::remove_unlikely_children(haplotypeStateNode* n, 
            double threshold) {
  for(size_t i = 0; i < n->number_of_children(); i++) {
    if(n->get_child(i)->prefix_likelihood() < threshold) {
      n->remove_child(n->get_child(i));
    }
  }
}

// TODO test []
haplotypeStateNode* haplotypeStateTree::alleles_to_state(
            const vector<alleleValue>& identifiers) const {
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

// TODO test []
vector<alleleValue> haplotypeStateTree::state_to_alleles(
            const haplotypeStateNode* state_node) const {
  const haplotypeStateNode* current_node = state_node;
  vector<alleleValue> backwards_allele_vector;
  while(current_node != root) {
    backwards_allele_vector.push_back(current_node->get_allele());
    current_node = current_node->get_parent();
  }
  vector<alleleValue> to_return;
  size_t upper = backwards_allele_vector.size() - 1;
  for(size_t i = 0; i < backwards_allele_vector.size(); i++) {
    to_return.push_back(backwards_allele_vector[upper - i]);
  }
  return to_return;
}

void haplotypeStateTree::set_initial_position(size_t position) {
  initial_position = position;
}
