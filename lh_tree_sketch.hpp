haplotypeStateNode{
  // the first allele which differentiates the node from the other children of its parents
  alleleValue local_identifier;
  vector<haplotypeStateNode*> children;
  haplotypeStateNode* parent;
  haplotypeMatrix* state;

  double partial_likelihood();
  void extend(alleleValue allele);
  void clear_state();

  haplotypeStateNode* create_child(alleleValue allele);
  void remove();

  vector<alleleValue> state_to_identifier();
};

// O(1) regardless of whether state->snapshot() has been called
double haplotypeStateNode::partial_likelihood() {
  return state->get_S();
}

// O(|minor allele here|) but not guranteed to produce a branch-able DP state until snapshot() is called
  void extend(alleleValue allele) {
}

haplotypeStateNode* haplotypeStateNode::create_child(alleleValue allele) {
  // this is O(|H|)
  state->snapshot();
  // copy snapshot of state to new state, add that as child
  state->extend(allele);
  // return pointer to child
}

void haplotypeStateNode::clear_state() {
  // delete the haplotypeMatrix* state of this node
  // to be called after all possible children of this node are built
}

void haplotypeStateNode::remove() {
  // if this is a leaf node, delete it and call remove() on its parent
}


vector<alleleValue> haplotypeStateNode::state_to_identifier() {
  // return sequence of local_identifiers of all parents up to root node
}

haplotypeStateNode* haplotypeStateTree::identifier_to_state(vector<alleleValue> identifier) {
  // use vector of alleleValues to traverse binary tree; return the final node reached
}
