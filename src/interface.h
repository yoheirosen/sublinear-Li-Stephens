#ifndef INTERFACE_H
#define INTERFACE_H

typedef long unsigned size_t;
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;

#ifdef __cplusplus
extern "C" {
#endif
  
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
  
haplotypeManager* haplotypeManager_build_abs_bound(
            char* reference_sequence,
            size_t ref_seq_length,
            size_t number_of_ref_sites,
            size_t* positions_of_ref_sites,
            size_t number_of_haplotypes,
            char* alleles_by_site_and_haplotype,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold);
            
haplotypeManager* haplotypeManager_build_interval_bound(
            char* reference_sequence,
            size_t ref_seq_length,
            size_t number_of_ref_sites,
            size_t* positions_of_ref_sites,
            size_t number_of_haplotypes,
            char* alleles_by_site_and_haplotype,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold);

// takes in an array of haplotypeStateNode*s of size 5. Indexed A-C-T-G-gap. 
// to this, it writes
//      pointers to the haplotypeStateNode*s corresponding to the children of
//      the state n, given that their prefix likelihood exceeds the threshold
//      for pruning
// We do not actually need to input the haplotypeManager itself because it is
// implied by the haplotypeStateNode *n being contained within it            
void haplotypeStateNode_get_next_options(
            haplotypeStateNode* n, 
            haplotypeStateNode** option_array);

// returns the conditional prefix likelihood for a state given the state of its
// prefix which contains one fewer site 
// Fast O(1) query given slow haplotypeManager pre-construction            
double haplotypeStateNode_local_probability(
            haplotypeStateNode* n, 
            haplotypeManager* hap_manager);
            
double haplotypeStateNode_total_probability(haplotypeStateNode* n);

// gets the allele of the haplotypeStateNode. Fast O(1) query
char haplotypeStateNode_allele(haplotypeStateNode* n);

// clears the haplotypeManager from memory
void haplotypeManager_delete(haplotypeManager* to_delete);

// this is the node which preceds the first "site." It is either empty, if the
// first "site" is also the first position within the subinterval. Otherwise it
// contains the prefix likelihood of the invariant interval which precedes it
haplotypeStateNode* haplotypeManager_get_root_node(
            haplotypeManager* hap_manager);

// step a state to the left by one "site"
haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n);

size_t haplotypeStateNode_number_of_children(haplotypeStateNode* n);

void haplotypeManager_print_transition_likelihoods(haplotypeManager* hap_manager);

void haplotypeManager_print_prefix_likelihoods(haplotypeManager* hap_manager);

#ifdef __cplusplus
}
#endif
#endif
