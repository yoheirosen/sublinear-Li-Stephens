#ifndef HAPLOTYPE_STATE_TREE_H
#define HAPLOTYPE_STATE_TREE_H

#include "probability.hpp"
#include "delay_multiplier.hpp"
#include "haplotype_state_node.hpp"
#include <unordered_set>
#include <cstdint>

using namespace std;

struct haplotypeStateTree{
private:
  siteIndex* reference;
  const penaltySet* penalties;
  const haplotypeCohort* cohort;

  // reference position at which the first tree-node starts
  size_t initial_position = SIZE_MAX;
  
public:
  haplotypeStateTree();
  haplotypeStateTree(
              siteIndex* reference, 
              const penaltySet* penalties, 
              const haplotypeCohort* cohort);
  ~haplotypeStateTree();
  
  haplotypeStateNode* root;
  
  void remove_node(haplotypeStateNode*& n);
  void remove_node_and_unshared_ancestors(haplotypeStateNode*& n);
  // void remove_unlikely_children(haplotypeStateNode* n, double threshold);
    
  haplotypeStateNode* alleles_to_state(
              const vector<alleleValue>& identifiers) const;
  vector<alleleValue> state_to_alleles(
              const haplotypeStateNode* state_node) const;

  void set_initial_position(size_t position);
};

#endif
