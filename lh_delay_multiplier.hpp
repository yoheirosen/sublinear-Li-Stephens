#ifndef DELAY_MULTIPLIER_H
#define DELAY_MULTIPLIER_H

#include <vector>
#include "lh_DP_map.hpp"

using namespace std;

// Shorthand for statements of complexity:
// |H|      number of haplotypes in population cohort
// n        length of haplotype since we last took a snapshot and reset the 
//          contents of the delayMap struct
// |slots|  a quantity whose expected value is a function of the frequency of
//          rare alleles and which is bounded by |H|
// M_avg    average across all sites of the number of haplotypes containing the 

// A delayMap is an O(|H| + n)-sized data structure which allows us to
// defer partial-probability update calculation of rows until the next site at
// which the row contains an allele seen in less than half the population
// This allows us to reduce the time complexity of the probability-calculation 
// DP to O(M_avg * n) from O(|H| * n) at the expense of a memory use increase to
// O(|H| + n) from O(|H|). However, the delayMap struct need not be
// stored and may be replaced with O(|H|) information as long as we call
// hard_update_all() first at at cost of O(|slots| + n) time
struct delayMap{
private:
  bool added_span;
  bool updated_maps;
  size_t dM_start;
  size_t current_site;
  size_t newest_index;
  vector<DPUpdateMap> maps_by_site;
  // This vector has size |H| and stores which map index the row
  // corresponds to
  vector<size_t> slots_by_row;
  // The following vectors all have equal |slots|
  vector<size_t> updated_to;
  vector<DPUpdateMap> maps_by_slot;
  vector<size_t> counts;
  // stores which map slots have been emptied so that new map
  // slots can be added in their place
  vector<size_t> empty_map_slots;
public:
  // steps forward the "current site" position
  void increment_site_marker();
  
  vector<size_t> rows_to_slots(const vector<size_t>& rows);
  
  // takes in a set of slot indices and extends their maps_by_slot
  // time complexity is O(|indices| + n)
  void update_maps(const vector<size_t>& indices);
  
  // Updates all maps to current position in preparation for taking a
  // "snapshot" of the current state of the DP.
  // This has a cost of O(|slots| + n) time but allows 
  // storage of the DP state in O(|population|) rather than 
  // O(|population| + |sites|) space
  void hard_update_all();
  
  void hard_clear_all();
  
  // un-associates row from slot and assigns an out-of-bounds slot index
  // and decrements slot
  void remove_row_from_slot(size_t row);
  // decrements slot row-count, destroys slot if row-count hits 0
  void decrement_slot(size_t slot);
  // clears slot and returns it to the list of empty slots
  void delete_slot(size_t slot);
  void add_identity_map();
  void add_map(double coefficient, double constant);
  void add_map(DPUpdateMap map);
  void assign_row_to_newest_index(size_t row);
  vector<size_t> get_map_indices();
  
  double get_coefficient(size_t row);  
  double get_constant(size_t row);
  DPUpdateMap get_map(size_t row);
  
  vector<DPUpdateMap> get_maps();
  
  vector<DPUpdateMap> get_maps_by_site();

  
  void update_map_with_span(DPUpdateMap span_map);
  void update_map_with_span(double coefficient, double constant);
  void add_map_for_site(DPUpdateMap site_map);
  void add_map_for_site(double coefficient, double constant);
  
  size_t last_update(size_t row);

  delayMap(size_t rows, size_t start = 0);
};

#endif