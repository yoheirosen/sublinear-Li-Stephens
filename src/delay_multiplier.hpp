#ifndef DELAY_MULTIPLIER_H
#define DELAY_MULTIPLIER_H

#include <vector>
#include "DP_map.hpp"
#include "row_set.hpp"

using namespace std;

typedef size_t eqclass_t;
typedef size_t row_t;
typedef size_t step_t;

struct mapHistory{
private:
	size_t start;
	vector<DPUpdateMap> elements;
  vector<size_t> previous;
  vector<DPUpdateMap> suffixes;
public:
  mapHistory();
  mapHistory(const DPUpdateMap& map, size_t start = 0);
  mapHistory(const mapHistory& other); 
  mapHistory(const mapHistory& other, size_t new_start);
  
  void reserve_length(size_t length);
	
	void push_back(const DPUpdateMap& map);
	
	DPUpdateMap& operator[](size_t i);
	DPUpdateMap& back();
  DPUpdateMap& suffix(size_t i);
  size_t& prev_site(size_t i);
  size_t start_site() const;
  
	size_t size() const;
  const vector<DPUpdateMap>& get_elements() const;
};

// Shorthand for statements of complexity:
// |H|      number of haplotypes in population cohort
// n        length of haplotype since we last took a snapshot and reset the 
//          contents of the lazyEvalMap struct
// |eqclass|  a quantity whose expected value is a function of the frequency of
//          rare alleles and which is bounded by |H|
// M_avg    average across all sites of the number of haplotypes containing the 

// A lazyEvalMap is an O(|H| + n)-sized data structure which allows us to
// defer partial-probability update calculation of rows until the next site at
// which the row contains an allele seen in less than half the population
// This allows us to reduce the time complexity of the probability-calculation 
// DP to O(M_avg * n) from O(|H| * n) at the expense of a memory use increase to
// O(|H| + n) from O(|H|). However, the lazyEvalMap struct need not be
// stored and may be replaced with O(|H|) information as long as we call
// hard_update_all() first at at cost of O(|eqclasses| + n) time
//
// TODO this is currently O(n) to copy. Speed it up. Though for single query 
// case doesn't matter if H < n(H^(2/3)). Certainly can assume that H < MAC
struct lazyEvalMap{
private:  
  step_t current_site = 0;
  step_t current_step = 0;
  mapHistory map_history;

  vector<eqclass_t> row_to_eqclass;                         // size = # haplotypes

  eqclass_t newest_eqclass = 0;
  vector<DPUpdateMap> eqclass_to_map;                    // size = # eqclasses
  vector<size_t> eqclass_size;                           // size = # eqclasses
  vector<step_t> eqclass_last_updated;                // size = # eqclasses
  // stores which map eqclasses have been emptied so that new map
  // eqclasses can be added in their place
  vector<eqclass_t> empty_eqclass_indices;                  // size = # eqclasses
	
	// Assignes a row to the eqclass containing the last DPUpdateMap added
	
	// un-associates row from eqclass and assigns an out-of-bounds eqclass index
  // and decrements eqclass
  // decrements eqclass row-count, destroys eqclass if row-count hits 0
  void decrement_eqclass(eqclass_t eqclass);
  // clears eqclass and returns it to the list of empty eqclasses
  void delete_eqclass(eqclass_t eqclass);
  vector<DPUpdateMap> suffixes;
    
  vector<size_t> site_n_classes;
  vector<eqclass_t> site_class_list_above;
  vector<eqclass_t> site_class_list_below;
  vector<eqclass_t> rep_eqclass_of_site;
  
  inline bool is_singleton(eqclass_t eqclass) const;
  inline bool is_front(eqclass_t eqclass) const;
  void site_class_list_pop(step_t site);
  eqclass_t get_eqclass_below(eqclass_t eqclass);
  eqclass_t get_eqclass_above(eqclass_t eqclass);
  eqclass_t get_rep_eqclass(step_t site);
  void set_eqclass_below(eqclass_t eqclass, eqclass_t below);
  void set_eqclass_above(eqclass_t eqclass, eqclass_t above);
  void set_rep_eqclass(step_t site, eqclass_t rep);
  void clear_rep_eqclass(step_t site);
  void eqclass_unset_last_updated(eqclass_t eqclass);
  void site_class_list_remove(size_t eqclass);
public:
  lazyEvalMap();
  lazyEvalMap(size_t rows, size_t start = 0);
  lazyEvalMap(const lazyEvalMap& other);
  
  void reserve_length(size_t length);
	
  void update_eqclass(eqclass_t eqclass);
  void remove_from_site(eqclass_t eqclass);
  void add_to_current_site(eqclass_t eqclass);
  
	void assign_row_to_newest_eqclass(row_t row);
	void remove_row_from_eqclass(row_t row);
	
  double evaluate(row_t row, double value) const;

  vector<size_t> rows_to_eqclasses(const rowSet& rows) const;
    
	void stage_map_for_site(const DPUpdateMap& site_map);
  void stage_map_for_span(const DPUpdateMap& span_map);

  // takes in a set of eqclass indices and extends their eqclass_to_map
  // time complexity is O(|indices| + n)
  void update_maps(const vector<eqclass_t>& eqclasses);
	
  void update_active_rows(const rowSet& active_rows);
  
  // Updates all maps to current position in preparation for taking a
  // "snapshot" of the current state of the DP.
  // This has a cost of O(|eqclasses| + n) time but allows 
  // storage of the DP state in O(|population|) rather than 
  // O(|population| + |sites|) space
  void hard_update_all();  
  // resets all maps to identity, lying in a single eqclass
  void hard_clear_all();

	void reset_rows(const rowSet& rows);

  // Adds a new eqclass containing the given DPUpdateMap
  void add_eqclass(const DPUpdateMap& map);
  void add_identity_eqclass();
	
  // get a vector of indices-among-eqclasses of maps assigned to rows
  const vector<size_t>& 			get_map_indices() const;
  const DPUpdateMap& 					get_map(row_t row) const;
	const vector<DPUpdateMap>& 	get_map_history() const;
  vector<DPUpdateMap>& 				get_maps();
  const vector<DPUpdateMap>& 	get_maps() const;
	double                      get_coefficient(row_t row) const;  
	double                      get_constant(row_t row) const;
    
  size_t number_of_eqclasses() const;
	size_t get_current_site() const;
  
  void increment_site_marker();

  size_t row_updated_to(row_t row) const;
	size_t last_update(row_t row) const;
  size_t get_eqclass(row_t row) const;
  
  void 
  condense_history(step_t top, step_t bottom);
};

#endif