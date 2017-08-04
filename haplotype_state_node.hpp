#ifndef HAPLOTYPE_STATE_NODE_H
#define HAPLOTYPE_STATE_NODE_H

#include "lh_probability.hpp"
#include <unordered_set>

using namespace std;

struct alleleAtSite{
  size_t site_index;
  alleleValue allele;
  alleleAtSite(size_t site, alleleValue allele);
  alleleAtSite();
};

struct haplotypeStateNode{
private:
  // tracks the last reference position
  size_t tail_length = 0;
  size_t suppressed_sites = 0;
  alleleAtSite unique_identifying_allele;
  haplotypeStateNode* parent = nullptr;
  vector<haplotypeStateNode*> children;
    
  // quickly delete members without preserving ordering
public:
  haplotypeMatrix* state = nullptr;  

  ~haplotypeStateNode();
  haplotypeStateNode();
  haplotypeStateNode(alleleAtSite unique_identifying_allele);
  haplotypeStateNode(alleleAtSite unique_identifying_allele, 
            haplotypeStateNode* parent);
  
  bool is_leaf();
  bool is_abandoned_stem();
  
  haplotypeStateNode* add_child(alleleAtSite a_at_s);
  haplotypeStateNode* add_child_copying_state(alleleAtSite a_at_s);
  void set_parent(haplotypeStateNode* n);
  
  haplotypeStateNode* get_child(alleleAtSite a_at_s) const;
  haplotypeStateNode* get_child(alleleValue a) const;
  haplotypeStateNode* get_child(size_t index) const;
  vector<haplotypeStateNode*> get_unordered_children() const;
  vector<haplotypeStateNode*> get_ordered_children();
  size_t number_of_children() const;
  haplotypeStateNode* get_parent() const;
  size_t node_to_child_index(haplotypeStateNode* n);
  
  void extend_state_by_site(alleleAtSite a_at_s);
  void extend_state_by_all_alleles();
  void extend_by_alleles_over_threshold(double threshold);
  
  void remove_child(haplotypeStateNode* c);
  void remove_child(alleleValue a);
  void remove_child_from_childvector(size_t index);
  void remove_child_from_childvector(haplotypeStateNode* n);
  
  void clear_state();
  void copy_state_from_node(haplotypeStateNode* other);  
  void compress_state();
  
  double prefix_likelihood() const;
  
  size_t get_end_site_index() const;  
  alleleValue identifying_allele() const;
  size_t identifying_site_index() const;
};

bool operator< (const haplotypeStateNode& n1, const haplotypeStateNode& n2);

#endif