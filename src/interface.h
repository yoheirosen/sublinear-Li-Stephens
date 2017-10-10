#ifndef INTERFACE_H
#define INTERFACE_H

typedef long unsigned size_t;
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;

#ifdef __cplusplus

extern "C" haplotypeManager* haplotypeManager_build(
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
            
extern "C" void haplotypeManager_get_next_options(
            haplotypeStateNode* n, 
            haplotypeStateNode** option_array);
            
extern "C" double haplotypeStateNode_local_probability(
            haplotypeStateNode* n, 
            haplotypeManager* hap_manager);
            
extern "C" double haplotypeStateNode_total_probability(haplotypeStateNode* n);

extern "C" char haplotypeStateNode_allele(haplotypeStateNode* n);

extern "C" void haplotypeManager_delete(haplotypeManager* to_delete);

extern "C" haplotypeStateNode* haplotypeManager_get_root_node(
            haplotypeManager* hap_manager);
            
extern "C" haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n);

extern "C" size_t haplotypeManager_number_of_children(haplotypeStateNode* n);

extern "C" void print_five(double* test);

#endif
#endif