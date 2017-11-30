#include "haplotype_state_node.hpp"
#include <algorithm>

using namespace std;

haplotypeStateNode::haplotypeStateNode() {
  
}

haplotypeStateNode::haplotypeStateNode(alleleValue allele) :
            allele(allele) {
  
}

haplotypeStateNode::haplotypeStateNode(alleleValue allele, 
          haplotypeStateNode* parent) :
          allele(allele), parent(parent) {
  
}

haplotypeStateNode::~haplotypeStateNode() {
  delete state;
  for(size_t i = 0; i < children.size(); i++) {
    delete children[i];
  }
}

haplotypeStateNode::haplotypeStateNode(const haplotypeStateNode& other) {
  S = other.S;
  allele = other.allele;
  parent = other.parent;
  children = {};
  state = new fastFwdAlgState(*state);
  marked_for_deletion = false;
}

bool haplotypeStateNode::is_root() const {
  return parent == nullptr;
}

bool haplotypeStateNode::is_leaf() const {
  return children.size() == 0;
}


bool haplotypeStateNode::is_abandoned_stem() const {
  return is_leaf() && state == nullptr;
}

haplotypeStateNode* haplotypeStateNode::add_child(alleleValue a) {
  haplotypeStateNode* new_child = new haplotypeStateNode(a, this);
  children.push_back(new_child);
  return new_child;
}

haplotypeStateNode* haplotypeStateNode::add_child_copying_state(
            alleleValue a) {
  // this->state->take_snapshot();
  haplotypeStateNode* new_child = add_child(a);
  new_child->copy_state_from_node(this);
  return new_child;
}

haplotypeStateNode* haplotypeStateNode::add_child_transferring_state(
            alleleValue a) {
  // this->state->take_snapshot();
  haplotypeStateNode* new_child = add_child(a);
  new_child->state = state;
  release_state();
  return new_child;
}

void haplotypeStateNode::set_parent(haplotypeStateNode* n) {
  parent = n;
}

haplotypeStateNode* haplotypeStateNode::get_child(size_t index) const {
  return children[index];
}

haplotypeStateNode* haplotypeStateNode::get_child(alleleValue a) const {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->allele == a) {
      return children[i];
    }
  }
  return nullptr;
}

size_t haplotypeStateNode::node_to_child_index(const haplotypeStateNode* c) const {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i] == c) {
      return i;
    }
  }
}

const vector<haplotypeStateNode*>& haplotypeStateNode::get_unordered_children() const {
  return children;
}

vector<haplotypeStateNode*> haplotypeStateNode::get_ordered_children() {
  sort(children.begin(), children.end());
  return children;
}

size_t haplotypeStateNode::number_of_children() const {
  return children.size();
}

haplotypeStateNode* haplotypeStateNode::get_parent() const {
  return parent;
}

void haplotypeStateNode::remove_child(haplotypeStateNode*& c) {
  remove_child_from_childvector(node_to_child_index(c));
  delete c;
  return;
}

void haplotypeStateNode::remove_child(alleleValue a) {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->get_allele() == a) {
      remove_child(i);
      return;
    }
  }
}

void haplotypeStateNode::remove_child(size_t index) {
  delete children[index];
  remove_child_from_childvector(index);      
  return;
}

void haplotypeStateNode::remove_child_from_childvector(size_t index) {
  children[index] = children.back();
  children.pop_back();
}

void haplotypeStateNode::remove_child_from_childvector(
            haplotypeStateNode* c) {
  remove_child_from_childvector(node_to_child_index(c));      
  return;
}

void haplotypeStateNode::clear_state() {
  if(state != nullptr) {
    S = state->prefix_likelihood();
  }
  delete state;
  state = nullptr;
}

void haplotypeStateNode::release_state() {
  if(state != nullptr) {
    S = state->prefix_likelihood();
  }
  state = nullptr;
}

// void haplotypeStateNode::delete_children_except(vector<size_t>& indices) {
//   vector<haplotypeStateNode*> temporary_children;
//   for(size_t i = 0; i < indices.size(); i++) {
//     temporary_children.push_back(children[i]);
//   }
// }


void haplotypeStateNode::copy_state_from_node(const haplotypeStateNode* other) {
  clear_state();
  state = new fastFwdAlgState(*(other->state));
}

void haplotypeStateNode::compress_state() {
  state->take_snapshot();
}

double haplotypeStateNode::prefix_likelihood() const {
  if(state != nullptr) {
    return state->prefix_likelihood();
  } else {
    return S;
  }
}

alleleValue haplotypeStateNode::get_allele() const {
  return allele;
}

double haplotypeStateNode::max_prefix_likelihood(const penaltySet* penalties) const {
  if(is_root()) {
    return 0;
  } else {
    return get_parent()->prefix_likelihood() - penalties->beta_value;
  }
}

void haplotypeStateNode::mark_for_deletion() {
  clear_state();
  marked_for_deletion = true;
}

bool haplotypeStateNode::is_marked_for_deletion() const {
  return marked_for_deletion;
}


bool operator< (const haplotypeStateNode& n1, const haplotypeStateNode& n2) {
  return(n1.prefix_likelihood() < n2.prefix_likelihood());
}