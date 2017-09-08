#ifndef SCORED_NODE_H
#define SCORED_NODE_H

#include "haplotype_state_node.hpp"
#include "probability.hpp"
#include "reference.hpp"

using namespace std;

// A concise representation of a haplotypeStateTree node; consists of pointer
// which allows access to the node's relationship as a sequence query to other
// sequence queries, as well as its probability DP matrix state given that it
// is being maintained
struct scoredNode{
private:
  const haplotypeStateNode* node;
  double local_probability = 0;
public:
  scoredNode();
  scoredNode(const haplotypeStateNode* node);
  // Should a child search state exist corresponding to the current node with 
  // sequence query extended by a single allele, the following methods return 
  // this node. If no such search state is found in the (trimmed) tree, these
  // return a null scoredNode 
  scoredNode extend_search(char c);
  scoredNode extend_search(alleleValue a);
  // Returns a node's parent. Equivalent to removing the terminal element from a
  // sequence query
  scoredNode step_back();
  
  double set_local_probability(const penaltySet* penalties);
  
  double get_score() const;
  double get_local_probability() const;
  alleleValue get_allele() const;
  
  size_t number_of_children() const;
  const haplotypeStateNode* get_node() const;
  scoredNode get_child(size_t i) const;
};

#endif