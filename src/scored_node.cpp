#include "scored_node.hpp"

scoredNode scoredNode::extend_search(char c) {
  return extend_search(char_to_allele(c));
}

scoredNode scoredNode::extend_search(alleleValue a) {
  haplotypeStateNode* result = node->get_child(a);
  return scoredNode(result);
}

scoredNode::scoredNode(const haplotypeStateNode* node) : node(node) {

}

scoredNode scoredNode::step_back() {
  return scoredNode(node->get_parent());
}

double scoredNode::get_score() const {
  return node->prefix_likelihood();
}

const haplotypeStateNode* scoredNode::get_node() const {
  return node;
}

alleleValue scoredNode::get_allele() const {
  return node->get_allele();
}

double scoredNode::get_local_probability() const {
  return local_probability;
}

double scoredNode::set_local_probability(const penaltySet* penalties) {
  if(!node->is_root()) {
    local_probability =
            node->prefix_likelihood() - node->max_prefix_likelihood(penalties);
  } else {
    local_probability = 0;
  }
}

size_t scoredNode::number_of_children() const {
  return node->number_of_children();
}

scoredNode scoredNode::get_child(size_t i) const {
  return scoredNode(node->get_child(i));
}