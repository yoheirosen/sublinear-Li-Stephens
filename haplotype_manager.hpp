#ifndef HAPLOTYPE_MANAGER_H
#define HAPLOTYPE_MANAGER_H

#include "lh_reference.hpp"
#include "haplotype_state_node.hpp"
#include "haplotype_state_tree.hpp"
#include <vector>
#include <string>

using namespace std;

struct haplotypeManager{
private:
  linearReferenceStructure* reference;
  haplotypeCohort* cohort;
  penaltySet* penalties;
  
  // reference-position of beginning of read
  size_t start_position;
  // reference-position of end of read
  size_t end_position;
  referenceSequence* reference_sequence;
  string* read_reference;
  
  // read-positions of {a_i} (all read-sites)
  vector<size_t> read_site_read_positions;
  
  vector<size_t> ref_site_below_read_site;
  void find_ref_sites_below_read_sites();
  
  // maps {a_i} -> {0,1}
  vector<bool> read_site_is_shared;
  // maps {c_k} -> {a_i}
  vector<size_t> shared_site_read_indices;
  void find_shared_sites();
  
  vector<size_t> subsequence_indices;
  void build_subsequence_indices();
    
  // sites and their read-reference alleles which occur between start_position
  // and the first shared site, if it exists, and end_position otherwise. This
  // too may be empty if there are zero reference sites
  vector<alleleAtSite> ref_sites_in_initial_span;
  // reference-sites occuring after each shared site; individual vectors may
  // be empty; vector-of-vectors will be empty if there are zero shared sites
  vector<vector<alleleAtSite> > ref_sites_after_shared_sites;
  void find_ref_only_sites_and_alleles();

  vector<size_t> invariant_penalties_by_read_site;
  void count_invariant_penalties();
  
  haplotypeStateTree tree;
  
  bool ref_sites = false;
  void check_for_ref_sites();
public:
  haplotypeManager(
          linearReferenceStructure* reference, haplotypeCohort* cohort, 
                penaltySet* penalties, referenceSequence* reference_sequence,
          vector<size_t> site_positions_within_read,
          string* read_reference = nullptr, 
                size_t start_reference_position = 0);
  
  size_t length();
  size_t read_sites();
  size_t shared_sites();
  
  size_t ref_position(size_t p);
  size_t read_position(size_t p);      
              
  size_t get_read_site_read_position(size_t i);
  size_t get_read_site_ref_position(size_t i);
  
  size_t index_among_shared_sites(size_t i);
  size_t index_among_read_only_sites(size_t i);
  
  size_t get_shared_site_read_index(size_t j);
  size_t get_shared_site_ref_index(size_t j);

  size_t get_ref_site_below_read_site(size_t i);
  
  float invariant_penalty_at_read_site(size_t i);
  
  bool contains_shared_sites();
  bool contains_ref_sites();
  bool contains_read_only_sites();
  
  // Tree functions
  
  // initializes haplotypeStateTree *tree; extends it to the positions before
  // the first shared site and returns the prefix likelihood up to this point
  void initialize_tree();
};


#endif