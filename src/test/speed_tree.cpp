#include <iostream>
#include <random>
#include <chrono>
#include "haplotype_manager.hpp"

int main(int argc, char* argv[]) {
  if(argc >= 2) {
    double cutoff = atof(argv[1]);
    size_t number_of_sites = 10;
    size_t number_of_haplotypes = 10;
    double alt_allele_frequency = 0.2;
    double shared_site_frequency = 0.2;
    if(argc >=3) {
      number_of_sites = strtoul(argv[2], NULL, 0);
    }
    if(argc >=4) {
      number_of_haplotypes = strtoul(argv[3], NULL, 0);
    }
    if(argc >= 6) {
      alt_allele_frequency = atof(argv[4]);
      shared_site_frequency = atof(argv[5]);
    }
        
    default_random_engine generator;
    generator.seed(chrono::system_clock::now().time_since_epoch().count());
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
    
    cout << "starting" << endl;
    
    for(size_t i = 0; i < number_of_sites; i++) {
      ref_sites.push_back(i);
      if(bernoulli_ref_is_shared(generator)) {
        shared_sites.push_back(i);
      } else {
        if(bernoulli_alt_allele(generator)) {
          size_t replacement_allele = which_allele(generator);
          read_alleles[i] = allele_to_char((alleleValue)replacement_allele);
        }
      }
    }
    
    cout << "ref\t\t" << ref_alleles << endl;
    cout << "read ref\t" << read_alleles << endl;
    string shared_locations = string(number_of_sites, ' ');
    for(size_t i = 0; i < shared_sites.size(); i++) {
      shared_locations[shared_sites[i]] = 'x';
    }
    cout << "shared sites\t" << shared_locations << endl << endl;
    
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
    
    for(size_t i = 0; i < number_of_haplotypes; i++) {
      if(i == 0) {
        cout << "haplotypes:\t";
      } else {
        cout << "\t\t";
      }
      for(size_t j = 0; j < number_of_sites; j++) {
        cout << allele_to_char(cohort_alleles[i][j]);
      }
      cout << endl;
    }
    
    if(shared_sites.size() == 0) {
      cout << "no tree to build" << endl;
      return 0;
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
    
    cout << endl << "building tree with " << hap_manager.shared_sites() << " shared sites ..." << endl;
    
    // hap_manager.set_cutoff_interval(cutoff);
    auto begin = chrono::high_resolution_clock::now();
    
    hap_manager.build_entire_tree_cutoff(cutoff);
      
    auto end = chrono::high_resolution_clock::now();
    
    cout << "built tree" << endl << endl;
    
    hap_manager.print_tree();
    // cout << endl;
    // hap_manager.print_tree_transitions();
    
    auto ms = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    cout << endl << "sites\t" << number_of_sites << "\tshared\t" << shared_sites.size() << "\t cutoff \t" << cutoff << "\tterminal leaves\t" << hap_manager.get_current_leaves().size() << "\t time \t" << ms << "\tms"<< endl;
    return 0;
  } else {
    cerr << "need cutoff for tree-pruning" << endl;
    return 1;
  }
}