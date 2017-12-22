#ifndef DELAY_MULTIPLIER_H
#define DELAY_MULTIPLIER_H

#include <vector>
#include "DP_map.hpp"
#include "row_set.hpp"

using namespace std;

struct mapHistory{
private:
	size_t start;
	vector<DPUpdateMap> elements;
public:
  mapHistory();
  mapHistory(const DPUpdateMap& map, size_t start = 0);
  mapHistory(const mapHistory& other); 
  mapHistory(const mapHistory& other, size_t new_start);
	
	void push_back(const DPUpdateMap& map);
	
	DPUpdateMap& operator[](size_t i);
	DPUpdateMap& back();
  
	size_t size() const;
  const vector<DPUpdateMap>& get_elements() const;
};

// Shorthand for statements of complexity:
// |H|      number of haplotypes in population cohort
// n        length of haplotype since we last took a snapshot and reset the 
//          contents of the lazyEvalMap struct
// |mapclass|  a quantity whose expected value is a function of the frequency of
//          rare alleles and which is bounded by |H|
// M_avg    average across all sites of the number of haplotypes containing the 

// A lazyEvalMap is an O(|H| + n)-sized data structure which allows us to
// defer partial-probability update calculation of rows until the next site at
// which the row contains an allele seen in less than half the population
// This allows us to reduce the time complexity of the probability-calculation 
// DP to O(M_avg * n) from O(|H| * n) at the expense of a memory use increase to
// O(|H| + n) from O(|H|). However, the lazyEvalMap struct need not be
// stored and may be replaced with O(|H|) information as long as we call
// hard_update_all() first at at cost of O(|mapclasses| + n) time
//
// TODO this is currently O(n) to copy. Speed it up. Though for single query 
// case doesn't matter if H < n(H^(2/3)). Certainly can assume that H < MAC
struct lazyEvalMap{
private:
	size_t current_site = 0;
  mapHistory maps_by_site;

  vector<size_t> row_to_mapclass;                         // size = # haplotypes

  size_t newest_mapclass = 0;
  vector<DPUpdateMap> mapclass_to_map;                    // size = # mapclasses
  vector<size_t> mapclass_size;                           // size = # mapclasses
  vector<size_t> mapclass_to_last_updated;                // size = # mapclasses
  // stores which map mapclasses have been emptied so that new map
  // mapclasses can be added in their place
  vector<size_t> empty_mapclass_indices;                  // size = # mapclasses
	
	// Assignes a row to the mapclass containing the last DPUpdateMap added
	
	// un-associates row from mapclass and assigns an out-of-bounds mapclass index
  // and decrements mapclass
  // decrements mapclass row-count, destroys mapclass if row-count hits 0
  void decrement_mapclass(size_t mapclass);
  // clears mapclass and returns it to the list of empty mapclasses
  void delete_mapclass(size_t mapclass);
public:
  lazyEvalMap();
  lazyEvalMap(size_t rows, size_t start = 0);
  lazyEvalMap(const lazyEvalMap& other);
	
	void assign_row_to_newest_mapclass(size_t row);
	void remove_row_from_mapclass(size_t row);
	
  double evaluate(size_t row, double value) const;

  vector<size_t> rows_to_mapclasses(const rowSet& rows) const;
	vector<bool> rows_to_mapclassmask(const rowSet& rows) const;
    
	void add_map_for_site(const DPUpdateMap& site_map);
  void update_map_with_span(const DPUpdateMap& span_map);

  // takes in a set of mapclass indices and extends their mapclass_to_map
  // time complexity is O(|indices| + n)
  void update_maps(const vector<size_t>& mapclasses);
	void update_maps(const vector<bool>& mapclassmask);
	
  void update_map_with_active_rows(const rowSet& active_rows);
  
  // Updates all maps to current position in preparation for taking a
  // "snapshot" of the current state of the DP.
  // This has a cost of O(|mapclasses| + n) time but allows 
  // storage of the DP state in O(|population|) rather than 
  // O(|population| + |sites|) space
  void hard_update_all();  
  // resets all maps to identity, lying in a single mapclass
  void hard_clear_all();

	void reset_rows(const rowSet& rows);

  // Adds a new mapclass containing the given DPUpdateMap
  void add_map(const DPUpdateMap& map);
  void add_identity_map();
	
  // get a vector of indices-among-mapclasses of maps assigned to rows
  const vector<size_t>& 			get_map_indices() const;
  const DPUpdateMap& 					get_map(size_t row) const;
	const vector<DPUpdateMap>& 	get_maps_by_site() const;
  vector<DPUpdateMap>& 				get_maps();
  const vector<DPUpdateMap>& 	get_maps() const;
	double                      get_coefficient(size_t row) const;  
	double                      get_constant(size_t row) const;
    
  size_t number_of_mapclasses() const;
	size_t get_current_site() const;
  
  void increment_site_marker();

  size_t row_updated_to(size_t row) const;
	size_t last_update(size_t row) const;
  size_t get_mapclass(size_t row) const;
	
	// statistics for benchmarking
	static size_t overall_traceback_steps;
	static size_t overall_slot_calculations;
	static size_t overall_vector_writes;
};

#endif