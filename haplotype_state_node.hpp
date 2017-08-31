#ifndef HAPLOTYPE_STATE_NODE_H
#define HAPLOTYPE_STATE_NODE_H

#include "lh_probability.hpp"
#include <unordered_set>

using namespace std;

struct haplotypeStateNode{
private:
  // tracks the last reference position
  alleleValue allele;
  haplotypeStateNode* parent = nullptr;
  vector<haplotypeStateNode*> children;
  
  double S;
  
public:
  // stores a probability calculation DP state
  haplotypeMatrix* state = nullptr;  

  ~haplotypeStateNode();
  haplotypeStateNode();
  haplotypeStateNode(alleleValue allele);
  haplotypeStateNode(alleleValue allele, haplotypeStateNode* parent);
  
  bool is_root() const;
  bool is_leaf() const;
  bool is_abandoned_stem() const;
  
  // adds a new child with allele marker a, and empty probability state
  haplotypeStateNode* add_child(alleleValue a);
  // adds a new child with probability state copied from this one; this state
  // is not yet extended to the new child's site
  haplotypeStateNode* add_child_copying_state(alleleValue a);

  void set_parent(haplotypeStateNode* n);
  
  haplotypeStateNode* get_child(alleleValue a) const;
  haplotypeStateNode* get_child(size_t index) const;
  size_t node_to_child_index(const haplotypeStateNode* n) const;

  vector<haplotypeStateNode*> get_unordered_children() const;
  vector<haplotypeStateNode*> get_ordered_children();

  size_t number_of_children() const;
  haplotypeStateNode* get_parent() const;
  
  void remove_child(haplotypeStateNode* c);
  void remove_child(alleleValue a);
  void remove_child_from_childvector(size_t index);
  void remove_child_from_childvector(haplotypeStateNode* n);
  
  void clear_state();
  void copy_state_from_node(const haplotypeStateNode* other);  
  void compress_state();
  
  double prefix_likelihood() const;
  double max_prefix_likelihood(penaltySet* penalties) const;
  
  alleleValue get_allele() const;
};

bool operator< (const haplotypeStateNode& n1, const haplotypeStateNode& n2);

#endif