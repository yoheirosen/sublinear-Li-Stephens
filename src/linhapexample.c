#include "interface.h"
#include <stdio.h>

int main() {
	char* reference_sequence = "AAAAAAA";
	long unsigned int ref_seq_length = 7;
	long unsigned int number_of_ref_sites = 3;
	long unsigned int positions_of_ref_sites[] = {0,2,5};
	long unsigned int number_of_haplotypes = 5;
	// all haplotypes are concatenated here... ie you have 5 haplotypes with 3 sites assigned alleles
	char alleles_by_site_and_haplotype[] = {"AAAAAAATAAATAGC"};
	double mutation_penalty = -9;
	double recombination_penalty = -6;
	long unsigned int read_DP_ref_start = 0;
	long unsigned int read_DP_site_count = 3;
	long unsigned int read_DP_site_offsets[] = {0,2,5};
	char* read_DP_sequence = "AAAAAAA";
	// prefix likelihood value below which to trim nodes of the tree
	// has to be < 0; value 0 denotes no threshold
	double threshold = 0;

	// builds and trims whole tree
	haplotypeManager* hap_manager = haplotypeManager_build_abs_bound(
		reference_sequence,
		ref_seq_length,
		number_of_ref_sites,
		positions_of_ref_sites,
		number_of_haplotypes,
		alleles_by_site_and_haplotype,
		mutation_penalty,
		recombination_penalty,
		read_DP_ref_start,
		read_DP_site_count,
		read_DP_site_offsets,
		read_DP_sequence,
		threshold);
	
	haplotypeStateNode* n = haplotypeManager_get_root_node(hap_manager);
	haplotypeStateNode* options[5];
	
	// fills options-vector with children of n; options vector must be
	// a minimum of number of children
	haplotypeStateNode_get_next_options(n, options);

	for(int i = 0; i < haplotypeStateNode_number_of_children(n); i++) {
		n = options[i];
		// what allele does this node have?		
		char allele = haplotypeStateNode_allele(n);
		printf("%c %lf\n", allele, haplotypeStateNode_local_probability(n, hap_manager));
	}
	
	// print the whole thing
	haplotypeManager_print_transition_likelihoods(hap_manager);
	printf("\n");
	haplotypeManager_print_prefix_likelihoods(hap_manager);
		
	return 0;
}
