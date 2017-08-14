#ifndef HAPLOTYPE_MANAGER_H
#define HAPLOTYPE_MANAGER_H

#include "lh_reference.hpp"
#include "haplotype_state_node.hpp"
#include "haplotype_state_tree.hpp"
#include "reference_sequence.hpp"
#include <vector>
#include <string>

using namespace std;

// A concise representation of a haplotypeStateTree node; consists of pointer
// which allows access to the node's relationship as a sequence query to other
// sequence queries, as well as its probability DP matrix state given that it
// is being maintained
struct scoredNode{
  double score;
  haplotypeStateNode* node;
  scoredNode(haplotypeStateNode* node, double score);
  scoredNode(haplotypeStateNode* node);
  // Should a child search state exist corresponding to the current node with 
  // sequence query extended by a single allele, the following methods return 
  // this node. If no such search state is found in the (trimmed) tree, these
  // return a null scoredNode 
  scoredNode extend_search(char c);
  scoredNode extend_search(alleleValue a);
  // Returns a node's parent. Equivalent to removing the terminal element from a
  // sequence query
  scoredNode pop_back();
};


// There are four indexing schemes:
// 1. Reference-position
// 2. Read-position
// 3. Reference-index
// 4. Read-index
// described as follows:
// 1.
// 2.
// 3.
// 4.
// We may convert between these by the following relationships

struct haplotypeManager{
private:
  linearReferenceStructure* reference;
  haplotypeCohort* cohort;
  penaltySet* penalties;
  
  // Reference-position of beginning of read
  size_t start_position;
  // Reference-position of end of read
  size_t end_position;
  
  // The full sequence of the reference, indexed by reference-position
  referenceSequence* reference_sequence;
  // The sequence which the reads consider to be "reference," indexed by offset
  // from the first position of the read-set
  string* read_reference;
  
  // read-positions of all read sites; maps {a_i} -> N
  vector<size_t> read_site_read_positions;
  
  // maps {a_i} -> {ref-sites}
  vector<size_t> ref_site_below_read_site;
  void find_ref_sites_below_read_sites();
  
  // maps {a_i} -> {0,1}; shared iff a_i \in {c_k}
  vector<bool> read_site_is_shared;
  // maps {c_k} -> {a_i}
  vector<size_t> shared_site_read_indices;
  void find_shared_sites();
  
  // maps {a_i} -> {b_j} \union {c_k}; sends {a_i}-indices to the respective
  // subsequence indices. This is possible because {b_j} and {c_k} are disjoint
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
  vector<size_t> invariant_penalties_by_ref_site;
  void count_invariant_penalties();
  
  haplotypeStateTree tree;
  
  // Stores whether there are reference sites spanned by our read-set. 
  // Check_for_ref_sites is O(log-length) so we want to avoid calling it
  bool ref_sites = false;
  void check_for_ref_sites();
  
  // Vector of pointers to nodes keep track of leaves of tree to be extended
  vector<haplotypeStateNode*> nodes_at_last_level_built;
  size_t last_level_built;
public:
  haplotypeManager(
          linearReferenceStructure* reference, haplotypeCohort* cohort, 
                penaltySet* penalties, referenceSequence* reference_sequence,
          vector<size_t> site_positions_within_read,
          string* read_reference, size_t start_reference_position);
  ~haplotypeManager();
  
  // Length in positions (ie base-pairs) of the region
  size_t length();
  // Number of read-sites from the read reference query
  size_t read_sites();
  // Number of sites which are both read-sites and reference-sites
  size_t shared_sites();
  
  // Converts read-position to ref_position
  size_t ref_position(size_t p);
  // Converts ref-position to read_position
  size_t read_position(size_t p);      
  
  // Convert from read site indices to positions            
  size_t get_read_site_read_position(size_t i);
  size_t get_read_site_ref_position(size_t i);
  
  
  size_t index_among_shared_sites(size_t i);
  size_t index_among_read_only_sites(size_t i);
  
  size_t get_shared_site_read_index(size_t j);
  size_t get_shared_site_ref_index(size_t j);

  size_t get_ref_site_below_read_site(size_t i);
  
  double invariant_penalty_at_read_site(size_t i);
  double invariant_penalty_at_ref_site(size_t i);
  
  bool contains_shared_sites();
  bool contains_ref_sites();
  bool contains_read_only_sites();
  
  // Tree functions
  
  // initializes haplotypeStateTree *tree; extends it to the positions before
  // the first shared site and returns the prefix likelihood up to this point
  void initialize_tree();
  
  void build_first_level(double threshold);
  void build_next_level(double threshold);
  void fill_in_level(double threshold, size_t start_site,
          size_t upper_bound_site); 
  void extend_final_level(double threshold);
  void build_entire_tree(double threshold);
  
  scoredNode find_node_by_prefix(string& prefix);
};


#endif