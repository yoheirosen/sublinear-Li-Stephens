#ifndef HAPLOTYPE_VCF_MANAGER_H
#define HAPLOTYPE_VCF_MANAGER_H

#include <vcflib>
#include "lh_reference.hpp"

using namespace std;

struct vcfManager{
  size_t num_phases;
  haplotypeCohort* cohort;
  linearReferenceStructure* reference;
  vcfManager(const char* vcf_path, const char* reference_sequence);
  vcfmanager(const string& vcf_path, const string& reference_sequence);
  void build_reference(const string& vcf_path, const string& reference_sequence);
};

#endif