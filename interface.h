#ifdef __cplusplus
#include "haplotype_manager.hpp"
#include "lh_reference.hpp"
#include "lh_probability.hpp"
#include "reference_sequence.hpp"
#include "vcf_manager.hpp"
#else
typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeStateNode haplotypeStateNode;
typedef struct penaltySet penaltySet;
#endif

using namespace std;

extern "C" haplotypeManager* haplotypeManager_build(
            char* reference_sequence,
            char* population_vcf_path,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold) {
  vcfManager vcf_manager(population_vcf_path, reference_sequence);
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         vcf_manager.num_phases);
  vector<size_t> read_sites(read_DP_site_offsets,
                            read_DP_site_offsets + read_DP_site_count);
  haplotypeManager* hap_manager = 
              new haplotypeManager(vcf_manager.reference, 
                                   vcf_manager.cohort, 
                                   penalties, 
                                   reference_sequence,
                                   read_sites,
                                   read_DP_sequence, 
                                   read_DP_ref_start);
  hap_manager->build_entire_tree(threshold);
  return hap_manager;
}

extern "C" void haplotypeManager_get_next_options(
            haplotypeManager* hap_manager,
            haplotypeStateNode* n, haplotypeStateNode** option_array) {
  size_t number_of_children = n->number_of_children();
  vector<haplotypeStateNode*> children = n->get_unordered_children();
  for(size_t i = 0; i < number_of_children; i++) {
    option_array[i] = children[i];
  }
  for(size_t i = number_of_children; i < 5; i++) {
    option_array[i] = nullptr;
  }
}

extern "C" double scoredNode_local_probability(haplotypeStateNode* n, penaltySet* penalties) {
  return n->prefix_likelihood() - n->max_prefix_likelihood(penalties);
}

extern "C" double haplotypeStateNode_total_probability(haplotypeStateNode* n) {
  return n->prefix_likelihood();
}

extern "C" char haplotypeStateNode_allele(haplotypeStateNode* n) {
  return allele_to_char(n->get_allele());
}

extern "C" void haplotypeManager_delete(haplotypeManager* to_delete) {
  delete to_delete->get_reference();
  delete to_delete->get_cohort();
  delete to_delete->get_penalties();
  delete to_delete;
}

extern "C" haplotypeStateNode* haplotypeManager_get_root_node(haplotypeManager* hap_manager) {
  return hap_manager->get_tree()->root;
}

extern "C" haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n) {
  return n->get_parent();
}