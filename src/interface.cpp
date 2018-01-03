#include "haplotype_manager.hpp"
#include "reference.hpp"
#include "probability.hpp"
#include "reference_sequence.hpp"
#include "input_haplotype.hpp"
#include "interface.h"
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

size_t haplotypeManager_get_num_shared_sites(haplotypeManager* hap_manager) {
  return hap_manager->shared_sites();
}

int haplotypeManager_read_index_is_shared(haplotypeManager* hap_manager, size_t read_site_index) {
  bool to_return = hap_manager->read_index_is_shared(read_site_index);
  return (int)to_return;
}

// double haplotypeManager_read_site_penalty(haplotypeManager* hap_manager, size_t read_site_index, char allele) {
//   if(hap_manager->read_matches(read_site_index, allele)) {
//     return penalties->one_minus_mu;
//   } else {
//     return penalties->mu;
//   }
// }

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
            double threshold) {
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         number_of_haplotypes);
  vector<alleleValue> ref_site_allele_vector;
  for(size_t i = 0; i < number_of_ref_sites; i++) {
    ref_site_allele_vector.push_back(allele::from_char(reference_sequence[positions_of_ref_sites[i]]));
  }
  vector<size_t> ref_site_position_vector = 
            vector<size_t>(positions_of_ref_sites, positions_of_ref_sites + number_of_ref_sites);
  
  siteIndex* reference = new siteIndex(ref_site_position_vector, ref_seq_length);  
    
  vector<vector<alleleValue> > haplotypes = 
            vector<vector<alleleValue> >(number_of_haplotypes,
                                         vector<alleleValue>(number_of_ref_sites, A));
                                         
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    for(size_t j = 0; j < number_of_ref_sites; j++) {
      haplotypes[i][j] = allele::from_char(alleles_by_site_and_haplotype[i*number_of_haplotypes + j]);
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
            double threshold) {
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         number_of_haplotypes);
  vector<alleleValue> ref_site_allele_vector;
  for(size_t i = 0; i < number_of_ref_sites; i++) {
    ref_site_allele_vector.push_back(allele::from_char(reference_sequence[positions_of_ref_sites[i]]));
  }
  vector<size_t> ref_site_position_vector = 
            vector<size_t>(positions_of_ref_sites, positions_of_ref_sites + number_of_ref_sites);
  
  siteIndex* reference =
            new siteIndex(ref_site_position_vector, 
                                         ref_seq_length);  
    
  vector<vector<alleleValue> > haplotypes = 
            vector<vector<alleleValue> >(number_of_haplotypes,
                                         vector<alleleValue>(number_of_ref_sites, A));
                                         
  for(size_t i = 0; i < number_of_haplotypes; i++) {
    for(size_t j = 0; j < number_of_ref_sites; j++) {
      haplotypes[i][j] = allele::from_char(alleles_by_site_and_haplotype[i*number_of_haplotypes + j]);
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
  hap_manager->build_entire_tree_interval(threshold);
  return hap_manager;
}

haplotypeManager* haplotypeManager_build_from_idx(
            char* reference_sequence,
            size_t ref_seq_length,
            siteIndex* reference,
            haplotypeCohort* cohort,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence) {
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         cohort->get_n_haplotypes());
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

  return hap_manager;
}

void haplotypeManager_build_tree_interval(haplotypeManager* hap_manager,
                                          double threshold) {
  hap_manager->build_entire_tree_interval(threshold);
}

void haplotypeStateNode_get_next_options(
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

haplotypeStateNode* haplotypeStateNode_get_child(
            haplotypeStateNode* parent,
            char allele) {
  return parent->get_child(allele::from_char(allele));
}

size_t haplotypeStateNode_number_of_children(haplotypeStateNode* n) {
  return n->number_of_children();
}

double haplotypeStateNode_local_probability(
            haplotypeStateNode* n, 
            haplotypeManager* hap_manager) {
  const penaltySet* penalties = hap_manager->get_penalties();
  return n->prefix_likelihood() - n->max_prefix_likelihood(penalties);
}

double haplotypeStateNode_total_probability(haplotypeStateNode* n) {
  return n->prefix_likelihood();
}

char haplotypeStateNode_allele(haplotypeStateNode* n) {
  return allele::to_char(n->get_allele());
}

void haplotypeManager_delete(haplotypeManager* to_delete) {
  delete to_delete->get_reference();
  delete to_delete->get_cohort();
  delete to_delete->get_penalties();
  delete to_delete;
}

haplotypeStateNode* haplotypeManager_get_root_node(
            haplotypeManager* hap_manager) {
  return hap_manager->get_tree()->root;
}

haplotypeStateNode* haplotypeStateNode_get_parent(haplotypeStateNode* n) {
  return n->get_parent();
}

void haplotypeManager_print_transition_likelihoods(haplotypeManager* hap_manager) {
  hap_manager->print_tree_transitions();
}

void haplotypeManager_print_prefix_likelihoods(haplotypeManager* hap_manager) {
  hap_manager->print_tree();
}

void haplotypeManager_print_terminal_nodes(haplotypeManager* hap_manager) {
  hap_manager->print_terminal_nodes();
}

haplotypeCohort* haplotypeCohort_init_empty(size_t number_of_haplotypes, siteIndex* ref) {
  return new haplotypeCohort(number_of_haplotypes, ref);
}

size_t haplotypeCohort_sum_MACs(haplotypeCohort* cohort) {
  return cohort->sum_MACs();
}

size_t haplotypeCohort_n_sites(haplotypeCohort* cohort) {
  return cohort->get_n_sites();
}

siteIndex* siteIndex_init_empty(size_t global_offset) {
  return new siteIndex(global_offset);
}

int64_t siteIndex_add_site(siteIndex* reference, size_t position) {
  return reference->add_site(position);
}

int haplotypeCohort_add_record(haplotypeCohort* cohort, size_t site) {
  return cohort->add_record(site);
}

int haplotypeCohort_set_sample_allele(
            haplotypeCohort* cohort, size_t site, size_t sample, char allele) {
  return cohort->set_sample_allele(site, sample, allele::from_char(allele));
}

void siteIndex_calc_spans(siteIndex* reference, size_t length) {
  reference->calculate_final_span_length(length);
}

void haplotypeCohort_populate_counts(haplotypeCohort* cohort) {
  cohort->populate_allele_counts();
}

void haplotypeCohort_delete(haplotypeCohort* cohort) {
  delete cohort;
}

void siteIndex_set_initial_span(siteIndex* ref, size_t length) {
  ref->set_initial_span(length);
}

size_t haplotypeCohort_n_haplotypes(haplotypeCohort* cohort) {
  return cohort->get_n_haplotypes();
}

size_t siteIndex_n_sites(siteIndex* reference) {
  return reference->number_of_sites();
}

int haplotypeManager_is_shared_site(haplotypeManager* hap_manager, size_t position) {
  return (int)(hap_manager->get_reference()->is_site(position));
}

void haplotypeManager_init_opt_idx(haplotypeManager* hap_manager,
                  								 char* r_alleles_1,
                  								 char* r_alleles_2) {
  char* ss_values_1 = (char*)malloc(hap_manager->shared_sites());
  char* ss_values_2 = (char*)malloc(hap_manager->shared_sites());
  for(size_t i = 0; i < hap_manager->shared_sites(); i++) {
    ss_values_1[i] = r_alleles_1[hap_manager->shared_index_to_read_index(i)];
    ss_values_2[i] = r_alleles_2[hap_manager->shared_index_to_read_index(i)];
  }
  hap_manager->set_option_index(ss_values_1, ss_values_2);
  free(ss_values_1);
  free(ss_values_2);
}

inputHaplotype* inputHaplotype_build(const char* ref_seq, 
                          const char* query, 
                          siteIndex* ref_struct,
                          size_t start_position) {
  inputHaplotype* to_return = new inputHaplotype(ref_seq, query, 
                                                 ref_struct, start_position, 
                                                 ref_struct->length_in_bp());
  return to_return;
}

void inputHaplotype_delete(inputHaplotype* in_hap) {
  delete in_hap;
}

fastFwdAlgState* fastFwdAlgState_initialize(siteIndex* reference,
                                            penaltySet* penalties,
                                            haplotypeCohort* cohort) {
  fastFwdAlgState* to_return = new fastFwdAlgState(reference, penalties, cohort);
  return to_return;
}

void fastFwdAlgState_delete(fastFwdAlgState* hap_matrix) {
  delete hap_matrix;
}

double fastFwdAlgState_score(fastFwdAlgState* hap_matrix, inputHaplotype* query) {
  double to_return = hap_matrix->calculate_probability(query);
  return to_return;
}

penaltySet* penaltySet_build(double recombination_penalty,
                             double mutation_penalty,
                             size_t number_of_haplotypes) {
  penaltySet* to_return = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         number_of_haplotypes);
  return to_return;
}

void penaltySet_delete(penaltySet* penalty_set) {
  delete penalty_set;  
}

typedef struct alleleVector alleleVector;
                              
slowFwdSolver* slowFwd_initialize(siteIndex* reference, penaltySet* penalties, haplotypeCohort* cohort) {
  return new slowFwdSolver(reference, penalties, cohort);
}

double slowFwd_solve_quadratic(slowFwdSolver* solver, inputHaplotype* q) {
  return solver->calculate_probability_quadratic(q->get_alleles(), q->get_start_index());
}

double slowFwd_solve_linear(slowFwdSolver* solver, inputHaplotype* q) {
  return solver->calculate_probability_linear(q->get_alleles(), q->get_start_index());
}

inputHaplotype* alleleVector_to_inputHaplotype(alleleVector* query, siteIndex* reference, size_t start_position, size_t end_position) {
  return new inputHaplotype(query->entries, vector<size_t>(query->size(), 0), reference, reference->start_position(), reference->length_in_bp());
}

inputHaplotype* haplotypeCohort_random_haplo(haplotypeCohort* cohort, siteIndex* reference, size_t generations, penaltySet* penalties, size_t length) {
  size_t start = reference->rand_interval_start(length);
  size_t start_site = reference->find_site_above(start);
  size_t end_site = reference->find_site_below(start + length);
  if(start_site <= end_site) {
    vector<alleleValue> random_haplo = cohort->rand_desc_haplo(generations, penalties->rho, penalties->mu, start, length);
    return new inputHaplotype(random_haplo, vector<size_t>(random_haplo.size(), 0), reference, start, length);
  } else {
    return new inputHaplotype();
  }
}

void alleleVector_delete(alleleVector* to_delete) {
  delete to_delete;
}

void slowFwdSolver_delete(slowFwdSolver* to_delete) {
  delete to_delete;
}

alleleVector* hC_separate_random(haplotypeCohort* cohort) {
  vector<size_t> choose_haplo = cohort->rand_haplos(1);
  alleleVector* to_return = new alleleVector(cohort->get_haplotype(choose_haplo[0]));
  haplotypeCohort* new_cohort = cohort->remove_haplotypes(choose_haplo);
  delete cohort;
  cohort = new_cohort;
  return to_return;
}

haplotypeCohort* hC_downsample_to(haplotypeCohort* cohort, size_t number) {
  return cohort->keep_haplotypes(number);
}

void aV_rebase_down(alleleVector* alleles, siteIndex* old_index, siteIndex* new_index) {
  alleleVector* new_alleles = new alleleVector;
  *new_alleles = allele::rebase_down(*alleles, *new_index);
  delete alleles;
  alleles = new_alleles;
  return;
}