#include "haplotype_manager.hpp"
#include "reference.hpp"
#include "probability.hpp"
#include "reference_sequence.hpp"
#include "interface.h"
#include <iostream>
#include <cmath>

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
      haplotypes[i][j] = char_to_allele(alleles_by_site_and_haplotype[i*number_of_haplotypes + j]);
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
      haplotypes[i][j] = char_to_allele(alleles_by_site_and_haplotype[i*number_of_haplotypes + j]);
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

haplotypeManager* haplotypeManager_build_int_from_index(
            char* reference_sequence,
            size_t ref_seq_length,
            linearReferenceStructure* reference,
            haplotypeCohort* cohort,
            double mutation_penalty, 
            double recombination_penalty,
            size_t read_DP_ref_start,
            size_t read_DP_site_count,
            size_t* read_DP_site_offsets,
            char* read_DP_sequence, 
            double threshold) {
  penaltySet* penalties = new penaltySet(recombination_penalty, 
                                         mutation_penalty, 
                                         cohort->get_haplotype_count());
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
  return allele_to_char(n->get_allele());
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

haplotypeCohort* haplotypeCohort_init_empty(size_t number_of_haplotypes, linearReferenceStructure* ref) {
  return new haplotypeCohort(number_of_haplotypes, ref);
}

linearReferenceStructure* linearReferenceStructure_init_empty(size_t global_offset) {
  return new linearReferenceStructure(global_offset);
}

int64_t linearReferenceStructure_add_site(
            linearReferenceStructure* reference, size_t position) {
  return reference->add_site(position);
}

int haplotypeCohort_add_record(haplotypeCohort* cohort, size_t site) {
  return cohort->add_record(site);
}

int haplotypeCohort_set_sample_allele(
            haplotypeCohort* cohort, size_t site, size_t sample, char allele) {
  return cohort->set_sample_allele(site, sample, char_to_allele(allele));
}

void linearReferenceStructure_calc_spans(linearReferenceStructure* reference, size_t length) {
  reference->calculate_final_span_length(length);
}

void haplotypeCohort_populate_counts(haplotypeCohort* cohort) {
  cohort->populate_allele_counts();
}

void linearReferenceStructure_set_initial_span(linearReferenceStructure* ref, size_t length) {
  ref->set_initial_span(length);
}


// Simulates an input from marginPhase consisting of:
//      1. positions of "sites"
//      2. consensus sequence at non-site positions
// In this particular simulation,
//      to generate 1., the parameter uncertainty_rate "U" is the probability 
//      that a site which is variable in the haplotype cohort is also variable, 
//      ie uncertain, from the perspective of the haplotype caller
//        - each reference cohort-site is chosen to be variable using a 
//          Bernoulli trial with p = U
//        - in order to have the same number of sites as the reference cohort,
//          additional sites are then chosen with uniform probability
// and
//       to generate the consensus sequence, a haplotype is copied, recombining
//       and mutating at random according to the probabilities from the 
//       penaltySet
// Return objects are passed by pointer and are:
//       return_read_sites -- the positions of sites "1"
//       n_return_read_sites -- number of sites
//       return_read_seq -- consensus sequence
// These can then be fed directly into the haplotypeManager constructor
void haplotypeCohort_sim_read_query(haplotypeCohort* cohort,
                                    const char* ref_seq,
                                    double mutation_rate,
                                    double recombination_rate,
                                    size_t cohort_size,
                                    double uncertainty_rate,
                                    size_t* return_read_sites,
                                    char* return_read_seq) {
  mutation_rate = pow(10, mutation_rate);
  recombination_rate = pow(10, recombination_rate);

  cohort->simulate_read_query(ref_seq,
                              mutation_rate,
                              recombination_rate,
                              uncertainty_rate,
                              return_read_sites,
                              return_read_seq);
}

size_t haplotypeCohort_n_haplotypes(haplotypeCohort* cohort) {
  return cohort->get_haplotype_count();
}

size_t linearReferenceStructure_n_sites(linearReferenceStructure* reference) {
  return reference->number_of_sites();
}