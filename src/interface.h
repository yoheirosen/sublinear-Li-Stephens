#ifndef INTERFACE_H
#define INTERFACE_H

typedef long unsigned size_t;
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;

#ifdef __cplusplus
extern "C" {
#endif
  
haplotypeManager* haplotypeManager_build(
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
            
void haplotypeStateNode_get_next_options(
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

size_t haplotypeStateNode_number_of_children(haplotypeStateNode* n);

void haplotypeManager_print(haplotypeManager* hap_manager);

#ifdef __cplusplus
}
#endif
#endif
