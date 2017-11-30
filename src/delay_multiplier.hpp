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
	DPUpdateMap& operator[](size_t i);
	DPUpdateMap& back();
	void push_back(const DPUpdateMap& map);
  size_t size() const;
  const vector<DPUpdateMap>& get_elements() const;
};

// Shorthand for statements of complexity:
// |H|      number of haplotypes in population cohort
// n        length of haplotype since we last took a snapshot and reset the 
//          contents of the lazyEvalMap struct
// |slots|  a quantity whose expected value is a function of the frequency of
//          rare alleles and which is bounded by |H|
// M_avg    average across all sites of the number of haplotypes containing the 

// A lazyEvalMap is an O(|H| + n)-sized data structure which allows us to
// defer partial-probability update calculation of rows until the next site at
// which the row contains an allele seen in less than half the population
// This allows us to reduce the time complexity of the probability-calculation 
// DP to O(M_avg * n) from O(|H| * n) at the expense of a memory use increase to
// O(|H| + n) from O(|H|). However, the lazyEvalMap struct need not be
// stored and may be replaced with O(|H|) information as long as we call
// hard_update_all() first at at cost of O(|slots| + n) time
struct lazyEvalMap{
private:
  bool added_span;
  bool updated_maps;
  size_t dM_start;
  size_t current_site = 0;
  size_t newest_index = 0;
  mapHistory maps_by_site;
  // This vector has size |H| and stores which map index the row
  // corresponds to
  vector<size_t> slots_by_row;
  // The following vectors all have equal |slots|
  vector<DPUpdateMap> maps_by_slot;
  vector<size_t> counts;
  vector<size_t> updated_to;
  // stores which map slots have been emptied so that new map
  // slots can be added in their place
  vector<size_t> empty_map_slots;
public:
  lazyEvalMap();
  lazyEvalMap(size_t rows, size_t start = 0);
  lazyEvalMap(const lazyEvalMap& other);
    
  // steps forward the "current site" position
  void increment_site_marker();
  
  vector<size_t> rows_to_slots(const vector<size_t>& rows) const;
  
  vector<size_t> rows_to_slots(const rowSet& rows) const;
  
  // takes in a set of slot indices and extends their maps_by_slot
  // time complexity is O(|indices| + n)
  void update_maps(const vector<size_t>& slots);
  
  // Updates all maps to current position in preparation for taking a
  // "snapshot" of the current state of the DP.
  // This has a cost of O(|slots| + n) time but allows 
  // storage of the DP state in O(|population|) rather than 
  // O(|population| + |sites|) space
  void hard_update_all();
  
  // resets all maps to identity, lying in a single slot
  void hard_clear_all();
  
  // un-associates row from slot and assigns an out-of-bounds slot index
  // and decrements slot
  void remove_row_from_slot(size_t row);
  // decrements slot row-count, destroys slot if row-count hits 0
  void decrement_slot(size_t slot);
  // clears slot and returns it to the list of empty slots
  void delete_slot(size_t slot);

  // Adds a new slot containing the given DPUpdateMap
  void add_map(double coefficient, double constant);
  void add_map(const DPUpdateMap& map);
  void add_identity_map();

  // Assignes a row to the slot containing the last DPUpdateMap added
  void assign_row_to_newest_index(size_t row);
  
  // get a vector of indices-among-slots of maps assigned to rows
  const vector<size_t>& get_map_indices() const;
  double get_coefficient(size_t row) const;  
  double get_constant(size_t row) const;
  const DPUpdateMap& get_map(size_t row) const;
  double evaluate(size_t row, double value) const;
  
  vector<DPUpdateMap>& get_maps();
  const vector<DPUpdateMap>& get_maps() const;
  
  const vector<DPUpdateMap>& get_maps_by_site() const;
  
  void update_map_with_span(const DPUpdateMap& span_map);
  void update_map_with_span(double coefficient, double constant);
  void add_map_for_site(const DPUpdateMap& site_map);
  void add_map_for_site(double coefficient, double constant);
  
  size_t last_update(size_t row);
  
  void reset_rows(const vector<size_t>& rows);
  void update_map_with_active_rows(const vector<size_t>& active_rows);
  
  void reset_rows(const rowSet& rows);
  void update_map_with_active_rows(const rowSet& active_rows);
    
  size_t number_of_slots() const;
  size_t row_updated_to(size_t row) const;
  size_t get_current_site() const;
  size_t get_slot(size_t row) const;
  
  void update_maps(const vector<bool>& slotmask);
  vector<bool> rows_to_slotmask(const rowSet& rows) const;
};

#endif