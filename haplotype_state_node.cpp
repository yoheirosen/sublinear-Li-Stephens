#include "haplotype_state_node.hpp"
#include <algorithm>

using namespace std;

alleleAtSite::alleleAtSite() {
  
}

alleleAtSite::alleleAtSite(size_t site, alleleValue allele) : site_index(site), 
          allele(allele) {
  
}

haplotypeStateNode::haplotypeStateNode() {
  
}

haplotypeStateNode::haplotypeStateNode(alleleAtSite unique_identifying_allele) :
            unique_identifying_allele(unique_identifying_allele) {
  
}

haplotypeStateNode::haplotypeStateNode(alleleAtSite unique_identifying_allele, 
          haplotypeStateNode* parent) :
          unique_identifying_allele(unique_identifying_allele), parent(parent) {
  
}

haplotypeStateNode::~haplotypeStateNode() {
  delete state;
  for(size_t i = 0; i < children.size(); i++) {
    delete children[i];
  }
}

bool haplotypeStateNode::is_leaf() {
  return (children.size() == 0);
}


bool haplotypeStateNode::is_abandoned_stem() {
  return (is_leaf() && state == nullptr);
}

haplotypeStateNode* haplotypeStateNode::add_child(alleleAtSite a_at_s) {
  haplotypeStateNode* new_child = new haplotypeStateNode(a_at_s);
  new_child->set_parent(this);
  children.push_back(new_child);
  return new_child;
}

haplotypeStateNode* haplotypeStateNode::add_child_copying_state(
            alleleAtSite a_at_s) {
  haplotypeStateNode* new_child = new haplotypeStateNode(a_at_s);
  new_child->set_parent(this);
  children.push_back(new_child);
  
  new_child->copy_state_from_node(this);
  new_child->extend_state_by_site(a_at_s);
  return new_child;
}


void haplotypeStateNode::set_parent(haplotypeStateNode* n) {
  parent = n;
}

haplotypeStateNode* haplotypeStateNode::get_child(size_t index) const {
  return children[index];
}

haplotypeStateNode* haplotypeStateNode::get_child(alleleAtSite a_at_s) const {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->unique_identifying_allele.allele == a_at_s.allele &&
              children[i]->unique_identifying_allele.site_index ==
              a_at_s.site_index) {
      return children[i];
    }
  }
  return nullptr;
}

haplotypeStateNode* haplotypeStateNode::get_child(alleleValue a) const {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->unique_identifying_allele.allele == a) {
      return children[i];
    }
  }
  return nullptr;
}

vector<haplotypeStateNode*> haplotypeStateNode::get_unordered_children() const {
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

size_t haplotypeStateNode::node_to_child_index(haplotypeStateNode* c) {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i] == c) {
      return i;
    }
  }
}

void haplotypeStateNode::extend_state_by_site(alleleAtSite a_at_s) {
  state->extend_probability_at_site(a_at_s.site_index, a_at_s.allele);
}

void haplotypeStateNode::extend_state_by_all_alleles(size_t site) {
  for(int i = 0; i < 5; i++) {
    alleleValue a = (alleleValue)i;
    haplotypeStateNode* new_child = 
              new haplotypeStateNode(alleleAtSite(site, a), this);
    new_child->copy_state_from_node(this);
    new_child->state->extend_probability_at_site(site, a);
  }
  delete this->state;
  this->state = nullptr;
}

void haplotypeStateNode::extend_by_alleles_over_threshold(size_t site,
            double threshold) {
  extend_state_by_all_alleles(site);
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->prefix_likelihood() < threshold) {
      delete children[i];
      remove_child_from_childvector(i);
      return;
    }
  }
}

void haplotypeStateNode::remove_child(haplotypeStateNode* c) {
  size_t i = node_to_child_index(c);
  delete children[i];
  remove_child_from_childvector(i);      
  return;
}

void haplotypeStateNode::remove_child(alleleValue a) {
  for(size_t i = 0; i < children.size(); i++) {
    if(children[i]->identifying_allele() == a) {
      delete children[i];
      remove_child_from_childvector(i);
      return;
    }
  }
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
  delete state;
}

void haplotypeStateNode::copy_state_from_node(haplotypeStateNode* other) {
  clear_state();
  state = new haplotypeMatrix(*(other->state));
}

void haplotypeStateNode::compress_state() {
  state->take_snapshot();
}

double haplotypeStateNode::prefix_likelihood() const {
  if(state != nullptr) {
    return state->prefix_likelihood();
  } else {
    return nan("");
  }
}

size_t haplotypeStateNode::get_end_site_index() const {
  return unique_identifying_allele.site_index + suppressed_sites;
}

alleleValue haplotypeStateNode::identifying_allele() const {
  return unique_identifying_allele.allele;
}

size_t haplotypeStateNode::identifying_site_index() const {
  return unique_identifying_allele.site_index;
}

bool operator< (const haplotypeStateNode& n1, const haplotypeStateNode& n2) {
  return(n1.prefix_likelihood() < n2.prefix_likelihood());
}