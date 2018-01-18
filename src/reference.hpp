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

//----------------------------------------------------------------------------------------------------------------------
// Index of site identifiers, positions, and the distances between them 
//----------------------------------------------------------------------------------------------------------------------
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
  int64_t add_site(size_t position);
  void calculate_final_span_length(size_t reference_length);

  //-- sizes -------------------------------------------------------------------
  size_t number_of_sites() const;
  size_t length_in_bp() const;

  //-- sites -------------------------------------------------------------------
  bool is_site(size_t actual_position) const;
  size_t get_site_index(size_t actual_position) const;
  size_t get_position(size_t site_index) const;
  
  //-- positions ---------------------------------------------------------------
  size_t start_position() const;
  size_t end_position() const;
  size_t pos_ref2global(size_t p) const;
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
  
  //-- random generators -------------------------------------------------------
  vector<alleleValue> make_child(const vector<alleleValue>& parent_0, const vector<alleleValue>& parent_1, double log_recomb_probability, double log_mutation_probability) const;
  vector<alleleValue> make_child(const vector<alleleValue>& parent_0, const vector<alleleValue>& parent_1, double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const;
  vector<size_t> rand_sites(size_t N) const;
  vector<size_t> rand_site_positions(size_t N) const;
  size_t rand_interval_start(size_t len) const;

  //-- downsampling ------------------------------------------------------------
  vector<alleleValue> downsample_to(const vector<alleleValue>& input, const siteIndex* index) const;
  
  //-- comparison --------------------------------------------------------------
  bool operator==(const siteIndex& other) const;
  bool operator!=(const siteIndex& other) const;
  bool operator<(const siteIndex& other) const;
  bool operator>(const siteIndex& other) const;
  bool operator<=(const siteIndex& other) const;
  bool operator>=(const siteIndex& other) const;
  vector<size_t> extra_sites(const siteIndex& other) const;
  
  //-- serialization -----------------------------------------------------------
  void write(std::ofstream& out) const;
  static siteIndex* read(std::ifstream& in);
};

//------------------------------------------------------------------------------
struct haplotypeCohort{
//------------------------------------------------------------------------------

private:
  typedef size_t haplo_id_t;
  typedef size_t site_idx_t;
  
  const siteIndex* reference;
  size_t number_of_haplotypes;
  bool finalized = false;

//------------------------------------------------------------------------------

  // TODO: ability to swap below for more succinct structures
  // maps [haplotypes] x [sites] -> alleles
  //      haplotype i           vector[i][ ]
  //      site j                vector[ ][j]
  vector<vector<alleleValue> > alleles_by_haplotype_and_site;

  // maps [sites] x [alleles] -> vectors of haplotype ids
  //      site i                vector[i][ ][ ]
  //      allele j              vector[ ][j][ ]
  //      haplotype rank k      vector[ ][ ][k]
  vector<vector<vector<haplo_id_t> > > haplotype_indices_by_site_and_allele;
  // TODO vector<vector<bv_tr_t> > > haplotype_indices_by_site_and_allele;
  vector<vector<rowSet> > active_rowSets_by_site_and_allele;

  // maps [sites] -> vectors of allele counts
  //      site i                vector[i][ ]
  //      allele j              vector[ ][j]
  vector<vector<size_t> > allele_counts_by_site_index;

//------------------------------------------------------------------------------
  haplotypeCohort* downsample_haplotypes(const vector<haplo_id_t>& ids, bool keep) const;

public:
//-- construction --

  haplotypeCohort(size_t cohort_size, const siteIndex* reference);
  haplotypeCohort(const vector<vector<alleleValue> >& haplotypes,
                  const siteIndex* reference);
  haplotypeCohort(const vector<string>& haplotypes, 
                  const siteIndex* reference);
  ~haplotypeCohort();
  
  // all-at-once
  void assign_alleles_at_site(site_idx_t i, vector<alleleValue> alleles_at_site);
  
  // per-site
  void set_column(const vector<alleleValue>& values);
  void set_column(const vector<alleleValue>& values, site_idx_t site);
  
  void add_record();
  void set_sample_allele(site_idx_t site, haplo_id_t sample, alleleValue a);
  
  void populate_allele_counts();
  rowSet build_active_rowSet(site_idx_t site, alleleValue a) const;
  
//-- basic attributes ----------------------------------------------------------

  const siteIndex* get_reference() const;  
  size_t get_n_sites() const;
  size_t get_n_haplotypes() const;
  
  bool operator==(const haplotypeCohort& other) const;
  bool operator!=(const haplotypeCohort& other) const;
  
//-- accessors -----------------------------------------------------------------
  
  // site x index -> allele
  alleleValue allele_at(site_idx_t site_index, haplo_id_t haplotype_index) const;
  const vector<alleleValue>& allele_vector_at_site(site_idx_t site_index) const;
  
  // index -> haplotype alleles
  const vector<alleleValue>& get_haplotype(haplo_id_t idx) const;

  // site x allele -> indices
  const vector<size_t>& get_matches(site_idx_t site_index, alleleValue a) const;
  vector<size_t> get_non_matches(site_idx_t site_index, alleleValue a) const;

  // site x allele -> counts
  bool match_is_rare(site_idx_t site_index, alleleValue a) const;
  size_t number_matching(site_idx_t site_index, alleleValue a) const;
  size_t number_not_matching(site_idx_t site_index, alleleValue a) const;

  // site -> mask
  vector<size_t> get_active_rows(site_idx_t site, alleleValue a) const;
  const rowSet& get_active_rowSet(site_idx_t site, alleleValue a) const;

//-- downsampling --------------------------------------------------------------

  // Removes sites where the minor allele frequency is below [double frequency].
  // For multiallelic sites, the minor allele frequency is taken to be the sum
  // of the minor allele frequencies with respect to the most common allele.
  // This constructs both a new (dynamically allocated) haplotypeCohort as well
  // as a new siteIndex, which is implicitly returned as
  haplotypeCohort* remove_sites_below_frequency(double frequency) const;

//-- more complex statistics ---------------------------------------------------

  alleleValue get_dominant_allele(site_idx_t site) const;
  size_t get_MAC(site_idx_t site) const;                                  // O(1)
  size_t sum_MACs() const;                                            // O(n)

//-- haplotype simulation ------------------------------------------------------
  // random generators
  vector<size_t> rand_haplos(size_t N) const;
  size_t rand_haplo_idx() const;
  size_t rand_haplo_idx(site_idx_t current) const;
  vector<alleleValue> rand_LS_haplo(double log_recomb_probability, double log_mutation_probability) const;
  vector<alleleValue> rand_LS_haplo(double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const;
  vector<alleleValue> rand_desc_haplo(size_t generations, double log_recomb_probability, double log_mutation_probability) const;
  vector<alleleValue> rand_desc_haplo(size_t generations, double log_recomb_probability, double log_mutation_probability, size_t start_site, size_t end_site) const;
  
//-- cohort editing ------------------------------------------------------------
  haplotypeCohort* remove_haplotypes(const vector<size_t>& ids_to_remove) const;
  haplotypeCohort* remove_haplotypes(size_t n_to_remove) const;
  haplotypeCohort* keep_haplotypes(const vector<size_t>& ids_to_keep) const;
  haplotypeCohort* keep_haplotypes(size_t n_to_keep) const;
  haplotypeCohort* remove_rare_sites(double max_rarity) const;
  
//-- serialization -------------------------------------------------------------
  void write(std::ofstream& out) const;
  void read_site(std::ifstream& in, site_idx_t site);
};

namespace haploRandom {
  vector<size_t> n_unique_uints(size_t N, size_t supremum);
  vector<size_t> n_unique_uints(size_t N, size_t supremum, const vector<size_t>* blacklist);
  alleleValue mutate(alleleValue a, double log_mutation_probability);
}

#endif