#include "scored_node.hpp"

scoredNode scoredNode::extend_search(char c) {
  return extend_search(char_to_allele(c));
}

scoredNode scoredNode::extend_search(alleleValue a) {
  haplotypeStateNode* result = node->get_child(a);
  return scoredNode(result);
}

scoredNode::scoredNode(haplotypeStateNode* node) : node(node) {
  if(node != nullptr) {
    if(node->state != nullptr) {
      *this = scoredNode(node, node->prefix_likelihood());
    } else {
      *this = scoredNode(node, 0);
    }
  } else {
    *this = scoredNode(nullptr, 0);
  }
}

scoredNode scoredNode::step_back() {
  return scoredNode(node->get_parent(), 0);
}


scoredNode::scoredNode(haplotypeStateNode* node, double score) : node(node),
            score(score) {
              
}

double scoredNode::get_score() const {
  return score;
}

haplotypeStateNode* scoredNode::get_node() const {
  return node;
}