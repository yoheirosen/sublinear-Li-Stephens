#ifndef DELAY_MULTIPLIER_H
#define DELAY_MULTIPLIER_H

#include <vector>
#include <stdexcept>
#include "DP_map.hpp"
#include "row_set.hpp"

using namespace std;

typedef size_t eqclass_t;
typedef size_t row_t;
typedef size_t step_t;

// TODO swap DPUpdateMap for pair of arrays
struct mapHistory{
private:
	vector<DPUpdateMap> maps;
  vector<step_t> previous;
  vector<step_t> next;
  // scratch area for suffix
  vector<DPUpdateMap> suffixes;
  vector<eqclass_t> rep_eqclasses;
  step_t leftmost;
public:
  class erased_error : public std::runtime_error {
    using std::runtime_error::runtime_error;
  };
  
  const static step_t CLEARED;
  const static step_t PAST_FIRST;
  const static eqclass_t NO_REP;
  const static step_t NOT_WRITEABLE;
  const static size_t NUMERICAL_UPPER_BOUND;
  
  mapHistory();
  // argument sets initial map
  mapHistory(const DPUpdateMap& map);
  void reserve(size_t length);
  
	size_t size() const;
  step_t rightmost_step() const;
  step_t leftmost_step() const;
  step_t past_end_step() const;
  
  // extend the mapHistory auxiliary vectors (previous, next) without constructing maps
  void increment_step();
  // construct maps. to be called after increment_step()
  void assign_rep_eqclass_to_new_step(eqclass_t eqclass);
  void assign_map_to_new_step(const DPUpdateMap& map);
  void postcompose_rightmost_map(const DPUpdateMap& map);
	
  const DPUpdateMap& operator[](step_t i) const;
	const DPUpdateMap& back() const;
  DPUpdateMap& suffix(step_t i);
  step_t previous_step(step_t i) const;
  step_t next_step(step_t i) const;
  eqclass_t rep_eqclass(step_t i) const;
  
  // bool step_cleared(step_t i) const;
  
  void set_rep_eqclass(step_t i, eqclass_t eqclass);
  void clear_rep_eqclass(step_t i);
  
  void set_new_leftmost(step_t i);
  void omit_step(step_t i);
  
  #ifdef DEBUG
    const vector<DPUpdateMap>& get_elements() const;
    vector<size_t> n_eqclasses;
    void dump_state(step_t marker) const;
    bool validate_step(step_t i) const;
    void print_index(size_t i) const;
  #endif
};

// Shorthand for statements of complexity:
// |H|      number of haplotypes in population cohort
// n        length of haplotype since we last took a snapshot and reset the 
//          contents of the delayedEvalMap struct
// |eqclass|  a quantity whose expected value is a function of the frequency of
//          rare alleles and which is bounded by |H|
// M_avg    average across all sites of the number of haplotypes containing the 

// A delayedEvalMap is an O(|H| + n)-sized data structure which allows us to
// defer partial-probability update calculation of rows until the next site at
// which the row contains an allele seen in less than half the population
// This allows us to reduce the time complexity of the probability-calculation 
// DP to O(M_avg * n) from O(|H| * n) at the expense of a memory use increase to
// O(|H| + n) from O(|H|). However, the delayedEvalMap struct need not be
// stored and may be replaced with O(|H|) information as long as we call
// hard_update_all() first at at cost of O(|eqclasses| + n) time

struct delayedEvalMap{
private:  
  step_t current_step = 0;
  
  // TODO implement
  vector<char> reused_active;
  vector<eqclass_t> reused_active_eqclasses;
  size_t n_active_eqclasses;
  
  mapHistory map_history;

  vector<eqclass_t> row_to_eqclass;                         // size = # haplotypes

  eqclass_t newest_eqclass = 0;
  vector<DPUpdateMap> eqclass_to_map;                    // size = # eqclasses
// TODO interleave
  vector<size_t> eqclass_size;                           // size = # eqclasses
  vector<step_t> eqclass_last_updated;                   // size = # eqclasses
  // stores which map eqclasses have been emptied so that new map
  // eqclasses can be added in their place
  vector<eqclass_t> empty_eqclass_indices;                  // size = # eqclasses
  
  const static eqclass_t INACTIVE_EQCLASS;
  const static size_t NUMERICAL_UPPER_BOUND;
  const static eqclass_t SINK_CLASS;
  
// TODO interleave
  vector<eqclass_t> eqclass_buddy_above;
  vector<eqclass_t> eqclass_buddy_below;
  
  bool step_has_singleton_eqclass(step_t i) const;
  bool eqclass_is_alone_at_step(eqclass_t eqclass) const;
#ifdef DEBUG
  bool validate_eqclass(eqclass_t eqclass) const;
  size_t count_eqclasses_at_step(step_t i) const;
  bool check_steps_valid(vector<step_t> steps) const;
  bool check_steps_valid() const;
#endif

  void push_back_at_nonempty_step(step_t i, eqclass_t eqclass);

  void delete_singleton_eqclass_from_step(eqclass_t eqclass, step_t i);
  void delete_non_singleton_eqclass_from_step(eqclass_t eqclass, step_t i);
  void delete_eqclass_from_step(eqclass_t eqclass);

  void move_singleton_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i);
  void move_non_singleton_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i);
  void move_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i);
  void move_leftmost_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i);

  void add_empty_eqclass_to_empty_step(const DPUpdateMap& map);

  void move_row_to_newest_eqclass(row_t row);
  void move_row_to_eqclass(row_t row, eqclass_t eqclass);
  void decrement_eqclass(eqclass_t eqclass);

  void update_eqclass(eqclass_t eqclass);
  void update_eqclasses(const vector<eqclass_t>& eqclasses);

public:
  delayedEvalMap();
  delayedEvalMap(size_t rows);
  
  void reserve_length(size_t length);
  
  void extend_by_new_step(const DPUpdateMap& map);
  void extend_value_only(const DPUpdateMap& map);
  
  double evaluate(row_t row, double value) const;
  void evaluate(const rowSet& active_rows, const vector<double>& values) const;
  void update(const rowSet& active_rows);
  
#ifdef TIME_PROBABILITY_INTERNALS
  void update_evaluate_and_move_rows(const rowSet& active_rows, vector<double>& values, double minority_correction, double* timer);
#else 
  void update_evaluate_and_move_rows(const rowSet& active_rows, vector<double>& values, double minority_correction);
#endif

  void hard_update_all();  
  void move_all_rows_to_newest_eqclass();

  // TODO efficient handling of small rowSets
  vector<eqclass_t> rows_to_eqclasses(const rowSet& rows) const;
  
  const DPUpdateMap& get_map(row_t row) const;
  
#ifdef DEBUG
  void print_eqclass(eqclass_t i) const;
  void dump_vectors() const;
  bool ensure_unique(eqclass_t eqclass, step_t allowed_step) const;
  bool ensure_deleted(eqclass_t eqclass) const;
#endif
};

#endif