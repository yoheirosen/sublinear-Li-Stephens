#include <iostream>
#include <random>
#include <chrono>
#include "reference.hpp"
#include "probability.hpp"
#include "input_haplotype.hpp"

using namespace std;

int main() {
  vector<double> allele_frequencies = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32};
  default_random_engine generator;
  poisson_distribution<int> poisson(33);
  for(int n = 0; n < allele_frequencies.size(); n++) {
    bernoulli_distribution bernoulli(allele_frequencies[n]);
    // replicates
    vector<size_t> lengths = {33, 66, 100, 333, 666, 1000, 3333, 6666, 10000};
    vector<size_t> populations = {5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000};
    for(int m = 0; m < 10; m++) {
      for(int l = 0; l < populations.size(); l++) {
        // lengths
        for(int k = 0; k < lengths.size(); k++) {
          // build reference
          size_t total_length = 0;
          size_t minor_allele_count = 0;
          vector<size_t> positions = {0};
          vector<alleleValue> reference_values (lengths[k], A);
          for(int j = 0; j < lengths[k] - 1; j++) {
            int next_interval = poisson(generator);
            positions.push_back(positions.back() + next_interval);
            total_length += next_interval;
          }
          linearReferenceStructure reference(positions, total_length, 
                      reference_values);
          vector<vector<alleleValue> > cohort_alleles;
          for(int i = 0; i < populations[l]; i++) {
            vector<alleleValue> haplotype_alleles;
            for(int j = 0; j < lengths[k]; j++) {
              if(bernoulli(generator)) {
                haplotype_alleles.push_back(T);
                minor_allele_count++;
              } else {
                haplotype_alleles.push_back(A);
              }
            }
            cohort_alleles.push_back(haplotype_alleles);
          }
          haplotypeCohort cohort(cohort_alleles, &reference);
          vector<size_t> augmentations (lengths[k] + 1, 0);
          inputHaplotype haplotype(reference_values, augmentations, &reference);
          penaltySet penalties(pow(10, -6), pow(10, -9), populations[l]);
          haplotypeMatrix hm(&reference, &penalties, &cohort);
          auto begin = chrono::high_resolution_clock::now();
          double prob = hm.calculate_probability(&haplotype);
          auto end = chrono::high_resolution_clock::now();
          auto ms = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
          cout << "frequency\t" << allele_frequencies[n] << "\tbp\t" << total_length << "\t sites \t" << lengths[k] << "\t population size \t" << populations[l] << "\t time \t" << ms << "\t ms\t" << minor_allele_count << "\t minor allele count" << endl;
        }
      }
    }
  }
  return 0;
}