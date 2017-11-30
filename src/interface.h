#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdint.h>

typedef long unsigned size_t;
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;
typedef struct linearReferenceStructure linearReferenceStructure;
typedef struct haplotypeCohort haplotypeCohort;
typedef struct inputHaplotype inputHaplotype;
typedef struct haplotypeMatrix haplotypeMatrix;

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
            
haplotypeManager* haplotypeManager_build_from_idx(
            char* reference_sequence,
            size_t ref_seq_length,
            linearReferenceStructure* reference,
            haplotypeCohort* cohort,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence);

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
            
haplotypeStateNode* haplotypeStateNode_get_child(
            haplotypeStateNode* parent,
            char allele);
            
double haplotypeStateNode_total_probability(haplotypeStateNode* n);

// gets the allele of the haplotypeStateNode. Fast O(1) query
char haplotypeStateNode_allele(haplotypeStateNode* n);

// clears the haplotypeManager from memory
void haplotypeManager_delete(haplotypeManager* to_delete);

int haplotypeManager_is_shared_site(haplotypeManager* hap_manager, size_t position);

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

void haplotypeManager_print_terminal_nodes(haplotypeManager* hap_manager);

size_t haplotypeManager_get_num_shared_sites(haplotypeManager* hap_manager);

int haplotypeManager_read_index_is_shared(haplotypeManager* hap_manager, size_t read_site_index);

// double haplotypeManager_read_site_penalty(haplotypeManager* hap_manager, size_t read_site_index, char allele);

// *****************************************************************************
// functions below used to construct haplotypeCohort and linearReferenceStructure
// from VCF
// *****************************************************************************

// initialize linearReferenceStructure with no sites
linearReferenceStructure* linearReferenceStructure_init_empty(size_t global_offset);

// initialize haplotypeCohort with no sites
haplotypeCohort* haplotypeCohort_init_empty(size_t number_of_haplotypes, linearReferenceStructure* ref);

// CONDITIONS:
//   must be called in ascending order of position, with no positions repeated
//   must be called before calling linearReferenceStructure_calc_spans, which 
//   locks the linearReferenceStructure
// Return value:
//   site index of just added site if conditions above are met 
//   -1 : positions added in non-ascending order
//   -2 : site with same position already
//   -3 : site-vector locked
int64_t linearReferenceStructure_add_site(
            linearReferenceStructure* reference, size_t position);

// Adds an empty record at the site specified, with all alleles set to unassigned
// CONDITIONS:
//   a site with the index specified must exist in the linearReferenceStructure
//   associated with the haplotypeCohort
//   this also must be the successor to the last site added to the haplotypeCohort
// Return value:
//   1 : successful
//   0 : locked
//   -1 : out of order
int haplotypeCohort_add_record(haplotypeCohort* cohort, size_t site);

// 1 : successful
// 0 : cohort locked
// -1 : already assigned
int haplotypeCohort_set_sample_allele(
            haplotypeCohort* cohort, size_t site, size_t sample, char allele);

// calculates final span length
// also locks the linearReferenceStructure from being modified
void linearReferenceStructure_calc_spans(linearReferenceStructure* reference, size_t length);

// builds allele-counts and lookup tables
// also locks the haplotypeCohort from being modified
void haplotypeCohort_populate_counts(haplotypeCohort* cohort);

void linearReferenceStructure_set_initial_span(linearReferenceStructure* ref, size_t length);

void haplotypeCohort_sim_read_query(haplotypeCohort* cohort,
                                    const char* ref_seq,
                                    double mutation_rate,
                                    double recombination_rate,
                                    size_t cohort_size,
                                    double uncertainty_rate,
                                    size_t* return_read_sites,
                                    char* return_read_seq);
                                    
void haplotypeCohort_sim_read_query_2(haplotypeCohort* cohort,
                                      const char* ref_seq,
                                      double mutation_rate,
                                      double recombination_rate,
                                      double uncertainty_rate,
                                      size_t** return_read_sites,
                                      size_t* n_read_sites,
                                      char* return_read_seq,
                                      char** r_s_alleles_1,
                                      char** r_s_alleles_2); 
                                    
size_t haplotypeCohort_n_haplotypes(haplotypeCohort* cohort);

size_t number_of_sites(haplotypeCohort* cohort);

size_t haplotypeCohort_sum_MACs(haplotypeCohort* cohort);

size_t linearReferenceStructure_n_sites(linearReferenceStructure* reference);

void haplotypeManager_init_opt_idx(haplotypeManager* hap_manager,
                  								 char* r_alleles_1,
                  								 char* r_alleles_2);

void haplotypeManager_build_tree_interval(haplotypeManager* hap_manager,
                                          double threshold);

inputHaplotype* inputHaplotype_build(const char* ref_seq, 
                          const char* query, 
                          linearReferenceStructure* ref_struct,
                          size_t start_position);

void inputHaplotype_delete(inputHaplotype* in_hap);

haplotypeMatrix* haplotypeMatrix_initialize(linearReferenceStructure* reference,
                                            penaltySet* penalties,
                                            haplotypeCohort* cohort);
                                              
void haplotypeMatrix_delete(haplotypeMatrix* hap_matrix);

double haplotypeMatrix_score(haplotypeMatrix* hap_matrix, inputHaplotype* query);

penaltySet* penaltySet_build(double recombination_penalty,
                             double mutation_penalty,
                             size_t number_of_haplotypes);

void penaltySet_delete(penaltySet* penalty_set);

#ifdef __cplusplus
}
#endif
#endif