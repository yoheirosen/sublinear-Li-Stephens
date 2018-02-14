#include <iostream>
#include <fstream>
#include "reference.hpp"

int main(int argc, char* argv[]) {
  if(argc != 2) {
    cerr << "serializer takes one argument, which is the vcf file path" << endl;
    return 1;
  }
  
  string vcf_path = argv[1];
  
  haplotypeCohort* temp = build_cohort(vcf_path);
  
  string slls_path = vcf_path.append(".slls");
  ofstream slls_out;
  slls_out.open(slls_path, ios::out | ios::trunc);
  temp->serialize_human(slls_out);
  slls_out.close();
  
  delete temp;
  
  return 0;
}