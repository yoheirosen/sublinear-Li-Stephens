#include <iostream>
#include <fstream>
#include <sys/time.h>
#include "reference.hpp"

int main(int argc, char* argv[]) {
  if(argc != 2) {
    cerr << "serializer takes one argument, which is the vcf file path" << endl;
    return 1;
  }

#ifdef TIME_COHORT_BUILD
  struct timeval build1, build2;
  gettimeofday(&build1, NULL);
#endif  
    
  string vcf_path = argv[1];
  haplotypeCohort* temp = build_cohort(vcf_path);

#ifdef TIME_COHORT_BUILD
  gettimeofday(&build2, NULL);
  double time_used_build = (double) (build2.tv_usec - build1.tv_usec) / 1000000 + (double) (build2.tv_sec - build1.tv_sec);
  cerr << vcf_path << "\t haplotypeCohort build time \t" << time_used_build << endl;
#endif

  string slls_path = vcf_path.append(".slls");
  ofstream slls_out;
  slls_out.open(slls_path, ios::out);

#ifdef TIME_COHORT_BUILD
  struct timeval write1, write2;
  gettimeofday(&write1, NULL);
#endif

  temp->serialize(slls_out);
  // temp->serialize_human(slls_out);
  slls_out.close();
  
#ifdef TIME_COHORT_BUILD
  gettimeofday(&write2, NULL);
  double time_used_output = (double) (write2.tv_usec - write1.tv_usec) / 1000000 + (double) (write2.tv_sec - write1.tv_sec);
  cerr << vcf_path << "\t haplotypeCohort file write time \t" << time_used_output << endl;
#endif
  
  delete temp;
  
  return 0;
}