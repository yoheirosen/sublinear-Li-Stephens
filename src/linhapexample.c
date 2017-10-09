#include "interface.h"

int main() {
	char* reference_sequence = "AAA";
	long unsigned int ref_seq_length = 3;
	long unsigned int number_of_ref_sites = 3;
	long unsigned int positions_of_ref_sites[] = {0,1,2};
	long unsigned int number_of_haplotypes = 2;
	char alleles_by_site_and_haplotype[] = {"ATACGC"};
	double mutation_penalty = 6;
	double recombination_penalty = 9;
	long unsigned int read_DP_ref_start = 0;
	long unsigned int read_DP_site_count = 3;
	long unsigned int read_DP_site_offsets[] = {0,1,2};
	char* read_DP_sequence = "CCC";
	double threshold = -30;

	haplotypeManager* hap_manager = haplotypeManager_build(
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
	haplotypeManager_get_next_options(n, options);
	n = options[0];
	haplotypeManager_get_next_options(n, options);
	n = options[1];
	double transition_probability =
		haplotypeStateNode_local_probability(n, hap_manager);
	char allele = haplotypeStateNode_allele(n);
	return 0;
}
