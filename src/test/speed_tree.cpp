#include <iostream>
#include <random>
#include <chrono>
#include "haplotype_manager.hpp"

int main() {
  double alt_allele_frequency = 0.;
  double shared_site_frequency = 0.1;
  size_t number_of_sites = 100;
  size_t number_of_haplotypes = 1000;
  
  double cutoff = -13;
  
  default_random_engine generator;
  bernoulli_distribution bernoulli_alt_allele(alt_allele_frequency);
  bernoulli_distribution bernoulli_ref_is_shared(shared_site_frequency);
  uniform_int_distribution<size_t> which_allele(1, 4);
  // poisson_distribution poisson_length(33);
  // size_t length_so_far = 0;
  // 
  vector<size_t> ref_sites;
  vector<size_t> shared_sites;
  string ref_alleles = string(number_of_sites, 'A');
  string read_alleles = string(number_of_sites, 'A');
  
  for(size_t i = 0; i < shared_sites.size(); i++) {
    ref_sites.push_back(i);
    if(bernoulli_ref_is_shared(generator)) {
      shared_sites.push_back(i);
    } else {
      if(bernoulli_alt_allele(generator)) {
        size_t replacement_allele = which_allele(generator);
        read_alleles[i] = (alleleValue)replacement_allele;
      }
    }
  }
  
  vector<vector<alleleValue> > cohort_alleles;
  for(int i = 0; i < number_of_haplotypes; i++) {
    vector<alleleValue> haplotype_alleles = vector<alleleValue>(number_of_sites, A);
    for(int j = 0; j < haplotype_alleles.size(); j++) {
      if(bernoulli_alt_allele(generator)) {
        size_t replacement_allele = which_allele(generator);
        haplotype_alleles[j] = (alleleValue)replacement_allele;
      }
    }
    cohort_alleles.push_back(haplotype_alleles);
  }
  
  linearReferenceStructure reference(ref_sites, ref_alleles);
  haplotypeCohort cohort(cohort_alleles, &reference);
  penaltySet penalties(-6, -9, number_of_haplotypes);
  haplotypeManager hap_manager = haplotypeManager(
    &reference, 
    &cohort, 
    &penalties, 
    ref_alleles.c_str(),
    shared_sites,
    read_alleles.c_str(), 
    0);

  auto begin = chrono::high_resolution_clock::now();
  
  hap_manager.build_entire_tree(cutoff);
    
  auto end = chrono::high_resolution_clock::now();
  
  // hap_manager.print();
  
  auto ms = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
  cout << "sites\t" << number_of_sites << "\t cutoff \t" << cutoff << "\t time \t" << ms << endl;
  return 0;
}