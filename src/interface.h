#ifdef __cplusplus
#include "haplotype_manager.hpp"
#include "reference.hpp"
#include "probability.hpp"
#include "reference_sequence.hpp"

using namespace std;
#else
typedef long unsigned size_t;
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;

haplotypeManager* haplotypeManager_build(
            char* reference_sequence,
            size_t number_of_ref_sites,
            size_t* positions_of_ref_sites,
            size_t number_of_haplotypes,
            char** alleles_by_site_and_haplotype,
            char* population_vcf_path,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold);
            
void haplotypeManager_get_next_options(
            haplotypeStateNode* n, 
            haplotypeStateNode** option_array);
            
double haplotypeStateNode_local_probability(
            haplotypeStateNode* n, 
            haplotypeManager* hap_manager);
            
double haplotypeStateNode_total_probability(haplotypeStateNode* n);

char haplotypeStateNode_allele(haplotypeStateNode* n);

void haplotypeManager_delete(haplotypeManager* to_delete);

haplotypeStateNode* haplotypeManager_get_root_node(
            haplotypeManager* hap_manager);
            
haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n);
#endif

#ifdef __cplusplus
// A haplotypeManager is the containter for all likelihood calculations
// it requires specification of:
//   - the base sequence of the reference
//   - penalties for mutations and recombinations (in log space)
//   - the data required to convert indexes between
//       - reference indexes
//       - and indexes relative to the start of the subinterval of the reference
//         which is being used
//              to build these, input of
//                  - read_DP_ref_start = the start of the subinterval under
//                    consideration in reference position 
//                  - read_DP_site_count = number of "sites" where we want to
//                    evaluate alleles choices for each haplotype
//                  - read_DP_site_offsets - the positions relative to the start
//                    of the subinterval
//   - the base sequence of the consensus sequence in the subinterval. This does
//     not, of course, need to be correctly specified at the "sites" because
//     that is precisely what we are trying to determine
//   - a cutoff threshold (in log-space)
// It also requires:
//   - a haplotype cohort. This can be read in directly from a vcf, to which the
//     path is provided
//       - the caveat in this case is that right now we can only handle SNPs or
//         single base deletions
//       - it might be preferable to take the intermediate step of reading the
//         vcf within marginPhase and then passing a 2D array of allele values
//         into this function instead
// From these, it builds out the prefix likelihoods at all possible sequences
// of values at "sites" which exceed the threshold value
//
// Performance is worst case exponential in the number of "sites" under
// consideration. This parameter must therefore be controlled with an
// appropriate choice of threshold. It otherwise depends linearly on:
//    - length of subinterval
//    - density of ref sites
//    - cohort size - average counts of non-major alleles across ref sites
extern "C" haplotypeManager* haplotypeManager_build(
            char* reference_sequence,
            size_t ref_seq_length,
            size_t number_of_ref_sites,
            size_t* positions_of_ref_sites,
            size_t number_of_haplotypes,
            char** alleles_by_site_and_haplotype,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold) {
  // TODO switch to direct read-in to cohort
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         number_of_haplotypes);
  vector<alleleValue> ref_site_allele_vector;
  for(size_t i = 0; i < number_of_ref_sites; i++) {
    ref_site_allele_vector.push_back(char_to_allele(reference_sequence[positions_of_ref_sites[i]]));
  }
  vector<size_t> ref_site_position_vector = 
            vector<size_t>(positions_of_ref_sites, positions_of_ref_sites + number_of_ref_sites);
  
  linearReferenceStructure* reference =
            new linearReferenceStructure(ref_site_position_vector, 
                                         ref_seq_length,
                                         ref_site_allele_vector);  
    
  vector<vector<alleleValue> > haplotypes = 
            vector<vector<alleleValue> >(number_of_haplotypes,
                                         vector<alleleValue>(number_of_ref_sites, A));
                                         
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    for(size_t j = 0; j < number_of_ref_sites; j++) {
      haplotypes[i][j] = char_to_allele(alleles_by_site_and_haplotype[i][j]);
    }
  }
  
  haplotypeCohort* cohort = 
            new haplotypeCohort(haplotypes, reference);  
  
  vector<size_t> read_sites(read_DP_site_offsets,
                            read_DP_site_offsets + read_DP_site_count);
  haplotypeManager* hap_manager = 
              new haplotypeManager(reference, 
                                   cohort, 
                                   penalties, 
                                   reference_sequence,
                                   read_sites,
                                   read_DP_sequence, 
                                   read_DP_ref_start);
  hap_manager->build_entire_tree(threshold);
  return hap_manager;
}

// takes in an array of haplotypeStateNode*s of size 5. Indexed A-C-T-G-gap. 
// to this, it writes
//      pointers to the haplotypeStateNode*s corresponding to the children of
//      the state n, given that their prefix likelihood exceeds the threshold
//      for pruning
// We do not actually need to input the haplotypeManager itself because it is
// implied by the haplotypeStateNode *n being contained within it
extern "C" void haplotypeManager_get_next_options(
            haplotypeStateNode* n, 
            haplotypeStateNode** option_array) {
  size_t number_of_children = n->number_of_children();
  vector<haplotypeStateNode*> children = n->get_unordered_children();
  for(size_t i = 0; i < number_of_children; i++) {
    option_array[i] = children[i];
  }
  for(size_t i = number_of_children; i < 5; i++) {
    option_array[i] = nullptr;
  }
}

// returns the conditional prefix likelihood for a state given the state of its
// prefix which contains one fewer site 
// Fast O(1) query given slow haplotypeManager pre-construction
extern "C" double haplotypeStateNode_local_probability(
            haplotypeStateNode* n, 
            haplotypeManager* hap_manager) {
  const penaltySet* penalties = hap_manager->get_penalties();
  return n->prefix_likelihood() - n->max_prefix_likelihood(penalties);
}

extern "C" double haplotypeStateNode_total_probability(haplotypeStateNode* n) {
  return n->prefix_likelihood();
}

// gets the allele of the haplotypeStateNode. Fast O(1) query
extern "C" char haplotypeStateNode_allele(haplotypeStateNode* n) {
  return allele_to_char(n->get_allele());
}

// clears the haplotypeManager from memory
extern "C" void haplotypeManager_delete(haplotypeManager* to_delete) {
  delete to_delete->get_reference();
  delete to_delete->get_cohort();
  delete to_delete->get_penalties();
  delete to_delete;
}

// this is the node which preceds the first "site." It is either empty, if the
// first "site" is also the first position within the subinterval. Otherwise it
// contains the prefix likelihood of the invariant interval which precedes it
extern "C" haplotypeStateNode* haplotypeManager_get_root_node(
            haplotypeManager* hap_manager) {
  return hap_manager->get_tree()->root;
}

// step a state to the left by one "site"
extern "C" haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n) {
  return n->get_parent();
}
#endif