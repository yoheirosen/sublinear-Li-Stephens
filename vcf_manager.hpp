#ifndef HAPLOTYPE_VCF_MANAGER_H
#define HAPLOTYPE_VCF_MANAGER_H

#include "Variant.h"
#include "reference.hpp"

using namespace std;

struct vcfManager{
  size_t num_phases;
  haplotypeCohort* cohort;
  linearReferenceStructure* reference;
  vcfManager(char* vcf_path, char* reference_sequence);
  vcfManager(string& vcf_path, string& reference_sequence);
  void build_reference(string& vcf_path, string& reference_sequence);
};

#endif