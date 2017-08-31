#include "haplotype_manager.hpp"
#include "lh_reference.hpp"
#include "lh_probability.hpp"
#include "reference_sequence.hpp"
#include "vcf_manager.hpp"

using namespace std;

typedef struct haplotypeManager haplotypeManager;
typedef struct haplotypeNode haplotypeNode;
typedef struct scoredNode scoredNode;

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
                                   reference_offset_of_read_DP);
  hap_manager->build_entire_tree(threshold);
  return hap_manager;
}

extern "C" scoredNode* haplotypeManager_get_next_options(
            haplotypeManager* hap_manager,
            scoredNode* n) {
  size_t number_of_children = n->number_of_children();
  vector<haplotypeStateNode*> children = n->get_unordered_children();
  scoredNode to_return[number_of_children];
  for(size_t i = 0; i < number_of_children; i++) {
    to_return[i] = scoredNode(children[i])
    if(local_probability == 0) {
      to_return[i].set_local_probability(penalties);
    }  
  }
  return to_return;
}

extern "C" double scoredNode_local_probability(scoredNode* n) {
  return n->get_local_probability();
}

extern "C" double scoredNode_total_probability(scoredNode* n) {
  return n->get_score();
}

extern "C" char scoredNode_allele(scoredNode* n) {
  return allele_to_char(n->get_allele());
}

extern "C" void haplotypeManager_delete(haplotypeManager* to_delete) {
  delete to_delete->reference;
  delete to_delete->cohort;
  delete to_delete->penalties;
  delete to_delete;
}

extern "C" scoredNode* haplotypeManager_get_root_node(haplotypeManager* hap_manager) {
  return scoredNode(hap_manager->get_tree()->root);
}

extern "C" scoredNode* scoredNode_get_parent(scoredNode* n) {
  return n->step_back();
}