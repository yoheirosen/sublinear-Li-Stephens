// A population reference cohort is expressed using two elements:
// 1. a siteIndex which stores positions of sites of variation
//    in a population
// 2. a haplotypeCohort which stores vectors of haplotype indices containing a
//    given allele at a given site. This is essential for rapid calculation of
//    the population likelihoods

#ifndef LINEAR_REFERENCE_STRUCTURE_H
#define LINEAR_REFERENCE_STRUCTURE_H

#include <string>
#include <vector>
#include <unordered_map>
#include "allele.hpp"
#include "row_set.hpp"

using namespace std;

//------------------------------------------------------------------------------
struct siteIndex{

private:
  size_t global_offset = 0;   // relative to the position on the chromosome
  size_t length = 0;          // length of spanned region in bp
  
  //-- site position data ------------------------------------------------------
  unordered_map<size_t, size_t> position_to_site_index;
  vector<size_t> site_index_to_position;
  
  //-- distances between sites or boundaries -----------------------------------
  vector<size_t> span_lengths;    // bp between sites, indexed by preceding site
                                  // final span is bp to end of region
  size_t leading_span_length;     // bp from beginning of region to first site

public:
  siteIndex(size_t global_offset);
  siteIndex(const vector<string>& haplotypes);
  siteIndex(const vector<size_t>& positions, size_t length);
  siteIndex(const vector<size_t>& positions, size_t length, size_t global_offset);
  ~siteIndex();
  
  //-- site-by-site construction -----------------------------------------------
  void set_initial_span(size_t length);
  // >= 0 : successful, returns site-index
  // -1 : out of order
  // -2 : collision
  int64_t add_site(size_t position);
  void calculate_final_span_length(size_t reference_length);

  //-- sizes -------------------------------------------------------------------
  size_t number_of_sites() const;
  size_t absolute_length() const;

  //-- sites -------------------------------------------------------------------
  bool is_site(size_t actual_position) const;
  size_t get_site_index(size_t actual_position) const;
  size_t get_position(size_t site_index) const;
  
  //-- positions ---------------------------------------------------------------
  size_t pos_ref2global(size_t p) const;
  size_t start_position() const;
  int64_t pos_global2ref(int64_t p) const;
  
  //-- spans -------------------------------------------------------------------
  bool has_span_before(size_t site_index) const;
  bool has_span_after(size_t site_index) const;
  size_t span_length_before(size_t site_index) const;
  size_t span_length_after(size_t site_index) const;
  
  //-- search ------------------------------------------------------------------
  // implemented as binary search
  // returns 1-past-the-end index if there is no site above
  size_t find_site_above(size_t position) const;
  // returns SIZE_MAX if there is no site below
  size_t find_site_below(size_t position) const;
};

//------------------------------------------------------------------------------
struct haplotypeCohort{
  
private:
  const siteIndex* reference;
  size_t number_of_haplotypes;
  bool finalized = false;

  //----------------------------------------------------------------------------
  // maps [haplotypes] x [sites] -> alleles
  //      haplotype i           vector[i][ ]
  //      site j                vector[ ][j]
  vector<vector<alleleValue> > haplotype_alleles_by_site_index;

  // maps [sites] x [alleles] -> vectors of haplotype ids
  //      site i                vector[i][ ][ ]
  //      allele j              vector[ ][j][ ]
  //      haplotype rank k      vector[ ][ ][k]
  vector<vector<vector<size_t> > > haplotype_indices_by_site_and_allele;

  // maps [sites] -> vectors of allele counts
  //      site i                vector[i][ ]
  //      allele j              vector[ ][j]
  vector<vector<size_t> > allele_counts_by_site_index;
public:
  haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
            const siteIndex* reference);
  haplotypeCohort(const vector<string>& haplotypes, 
            const siteIndex* reference);
  haplotypeCohort(size_t cohort_size, const siteIndex* reference);
  ~haplotypeCohort();
  
  void populate_allele_counts();
  
  // Removes sites where the minor allele frequency is below [double frequency].
  // For multiallelic sites, the minor allele frequency is taken to be the sum
  // of the minor allele frequencies with respect to the most common allele.
  // This constructs both a new (dynamically allocated) haplotypeCohort as well
  // as a new siteIndex, which is implicitly returned as

  haplotypeCohort* remove_sites_below_frequency(double frequency) const;
  const siteIndex* get_reference() const;  
  
  void assign_alleles_at_site(size_t i, vector<alleleValue> alleles_at_site);
  
  size_t size() const;
  size_t get_haplotype_count() const;
  alleleValue get_dominant_allele(size_t site) const;
  size_t get_MAC(size_t site) const;
  size_t sum_MACs() const;
  
  const vector<size_t>& get_matches(size_t site_index, alleleValue a) const;
  vector<size_t> get_non_matches(size_t site_index, alleleValue a) const;
  alleleValue allele_at(size_t site_index, size_t haplotype_index) const;
  size_t number_matching(size_t site_index, alleleValue a) const;
  size_t number_not_matching(size_t site_index, alleleValue a) const;
  bool match_is_rare(size_t site_index, alleleValue a) const;
  
  const vector<alleleValue>& get_alleles_at_site(size_t site_index) const;
  
  vector<size_t> get_active_rows(size_t site, alleleValue a) const;
  rowSet get_active_rowSet(size_t site, alleleValue a) const;

  // 1 : successful
  // 0 : locked
  // -1 : out of order
  int add_record(size_t site);
  // 1 : successful
  // 0 : cohort locked
  // -1 : already assigned
  int set_sample_allele(size_t site, size_t sample, alleleValue a);
  
  void set_column(const vector<alleleValue>& values);
  void set_column(const vector<alleleValue>& values, size_t site);
  
  void simulate_read_query(const char* ref_seq,
                           double mutation_rate,
                           double recombination_rate,
                           double uncertainty_rate,
                           size_t* return_read_sites,
                           char* return_read_seq) const;
                           
  void simulate_read_query_2(
                           const char* ref_seq,
                           double mutation_rate,
                           double recombination_rate,
                           double uncertainty_rate,
                           size_t** return_read_sites,
                           size_t* n_read_sites,
                           char* return_read_seq,
                           char** r_s_alleles_1,
                           char** r_s_alleles_2) const; 
};

#endif