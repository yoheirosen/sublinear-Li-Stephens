#include "vcf_manager.hpp"
#include <string>

vcfManager::vcfManager(const string& vcf_path, 
                       const string& reference_sequence) {
  build_reference(vcf_path, reference_sequence);
}

vcfManager::vcfManager(const char* vcf_path,
                       const char* reference_sequence) {
  const string path_string = string(vcf_path);
  const string ref_string = string(reference_sequence);
  build_reference(path_string, ref_string);
}

void vcfManager::build_reference(const string& vcf_path, 
                            const string& reference_sequence) {
  vcflib::VariantCallFile variant_file;
    
  variant_file.open(vcf_path);
  if(variant_file.is_open()) {
    size_t num_samples = variant_file.sampleNames.size();
    num_phases = num_samples * 2;
    
    vector<size_t> site_positions;
    vcflib::Variant var(variant_file);
    
    while(variant_file.getNextVariant(var)) {
      site_positions.push_back(var.position - 1);
    }
    
    reference = new linearReferenceStructure(site_positions, reference_sequence);

    var = Variant(variant_file);
    alleleValue allele;
    size_t current_site = 0;
    
    // to read in to the haplotypeCohort initializer
    vector<vector<alleleValue> > haplotypes = 
              vector<vector<alleleValue> >(num_phases,
                        vector<alleleValue>(site_positions.size(), (alleleValue)4));
    
    while(variant_file.getNextVariant(var)) {
      for(size_t sample = 0; sample < num_samples; sample++) {
        string& sample_name = variant_file.sampleNames[sample];
        string genotype = var.getGenotype(sample_name);
        auto bar_pos = genotype.find('|');
        vector<string> alt_indices({genotype.substr(0, bar_pos),
                genotype.substr(bar_pos + 1)});
        
        for(size_t phase_offset = 0; phase_offset < 2; phase_offset++) {
          string& alt_string = alt_indices[phase_offset];
          if (alt_string == ".") {
            allele = reference->get_reference_allele_at_site(current_site);
          } else {
            int all_value = stoi(alt_string);
            if(all_value == 0) {
              allele = reference->get_reference_allele_at_site(current_site);
            } else {
              allele = str_to_allele(var.alleles[all_value]);
            }
          }
          haplotypes[current_site][sample*2 + phase_offset] = allele;
        }
      }
      current_site++;
    }
    
    cohort = new haplotypeCohort(haplotypes, reference);
  }  
}