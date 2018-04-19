#include "row_set.hpp"
#include "delay_multiplier.hpp"
#include "math.hpp"

#include <iostream>

#ifdef DEBUG
#include <cassert>
#endif

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// -- mapHistory -------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

const step_t mapHistory::CLEARED = SIZE_MAX - 4;
const step_t mapHistory::PAST_FIRST = SIZE_MAX - 3;
const eqclass_t mapHistory::NO_REP = SIZE_MAX - 2;
const step_t mapHistory::NOT_WRITEABLE = SIZE_MAX - 5;
const size_t mapHistory::NUMERICAL_UPPER_BOUND = SIZE_MAX - 6;

const eqclass_t delayedEvalMap::INACTIVE_EQCLASS = SIZE_MAX - 1;
const size_t delayedEvalMap::NUMERICAL_UPPER_BOUND = SIZE_MAX - 6;
const eqclass_t delayedEvalMap::SINK_CLASS = 0;

mapHistory::mapHistory() {
}

mapHistory::mapHistory(const DPUpdateMap& map) {
	maps = {map};
  previous = {PAST_FIRST, 0};
  next = {1, 2};
  suffixes = {DPUpdateMap::IDENTITY};
  rep_eqclasses = {1};
  leftmost = 0;
  #ifdef DEBUG
    n_eqclasses = {1};
  #endif
}

void mapHistory::reserve(size_t length) {
  maps.reserve(length);
  suffixes.reserve(length);
  previous.reserve(length + 1);
  next.reserve(length + 1);
  rep_eqclasses.reserve(length);
}

step_t mapHistory::size() const {
  return rep_eqclasses.size();
}

step_t mapHistory::rightmost_step() const {
  return rep_eqclasses.size() - 1;
}

step_t mapHistory::leftmost_step() const {
  return leftmost;
}

step_t mapHistory::past_end_step() const {
  return rep_eqclasses.size();
}

void mapHistory::increment_step() {
#ifdef DEBUG
  if(maps.size() != previous.size() - 1) {
    cerr << "aux vector size " << previous.size() << "; map vector size " << maps.size() << endl;
    throw std::runtime_error("mapHistory::increment_step() : step incremented out of order");
  }
#endif
  previous.push_back(previous.size() - 1);
  next.push_back(next.size() + 1);
}

void mapHistory::assign_rep_eqclass_to_new_step(eqclass_t eqclass) {
#ifdef DEBUG
  assert(eqclass != 0);
#endif
  rep_eqclasses.push_back(eqclass);
}

void mapHistory::assign_map_to_new_step(const DPUpdateMap& map) {
  #ifdef DEBUG
    if(maps.size() != previous.size() - 2) {
      cerr << "aux vector size " << previous.size() << "; map vector size " << maps.size() << endl;
      throw std::runtime_error("mapHistory::assign_map_to_new_step(map) : step incremented out of order");
    }
  #endif
  suffixes.push_back(map);
  maps.push_back(map);
}

void mapHistory::postcompose_rightmost_map(const DPUpdateMap& map) {
  maps.back() = map.of(maps.back());
}

const DPUpdateMap& mapHistory::operator[](step_t i) const {
#ifdef DEBUG
  assert(validate_step(i));
#endif
  return maps[i];
}

const DPUpdateMap& mapHistory::back() const {
	return maps.back();
}

DPUpdateMap& mapHistory::suffix(step_t i) {
#ifdef DEBUG
  assert(validate_step(i));
#endif
  return suffixes[i];
}

step_t mapHistory::previous_step(step_t i) const {
#ifdef DEBUG
  assert(validate_step(i));  
#endif
  return previous[i];
}

step_t mapHistory::next_step(step_t i) const {
#ifdef DEBUG
  assert(validate_step(i));  
#endif
  return next[i];
}

eqclass_t mapHistory::rep_eqclass(step_t i) const {
#ifdef DEBUG
  assert(validate_step(i));
  assert(rep_eqclasses[i] != 0);
#endif
  return rep_eqclasses[i];
}

void mapHistory::set_rep_eqclass(step_t i, eqclass_t eqclass) {
#ifdef DEBUG
  assert(validate_step(i));
  assert(rep_eqclasses[i] != CLEARED);
#endif
  rep_eqclasses[i] = eqclass;
}

void mapHistory::clear_rep_eqclass(step_t i) {
#ifdef DEBUG
  assert(validate_step(i));
#endif
  rep_eqclasses[i] = NO_REP;
}

void mapHistory::set_new_leftmost(step_t i) {
#ifdef DEBUG
  assert(validate_step(i));
  previous[leftmost] = CLEARED;
  next[leftmost] = CLEARED;
#endif
  // this is the minimal work needed
  previous[i] = PAST_FIRST;
  leftmost = i;
}

void mapHistory::omit_step(step_t i) {
#ifdef DEBUG
  assert(i < rightmost_step());
  assert(validate_step(i));
  assert(i != leftmost);
  clear_rep_eqclass(i);
#endif
  step_t this_next = next[i];
  step_t this_previous = previous[i];
#ifdef DEBUG
  next[i] = CLEARED;
  previous[i] = CLEARED;
#endif
  next[this_previous] = this_next;
  previous[this_next] = this_previous;
  maps[this_previous] = maps[i].of(maps[this_previous]);
}

#ifdef DEBUG
  void mapHistory::print_index(size_t i) const {
    if(i == CLEARED) {
      cerr << "CLRD";
    } else if(i == PAST_FIRST) {
      cerr << "PST1";
    } else if(i == NO_REP) {
      cerr << "NREP";
    } else {
      cerr << i;
    }
  }
  
  void delayedEvalMap::print_eqclass(eqclass_t i) const {
    if(i == INACTIVE_EQCLASS) {
      cerr << "INAC";
    } else if(i == SINK_CLASS) {
      cerr << "SINK";
    } else {
      cerr << i;
    }
  }
  
  void delayedEvalMap::dump_vectors() const {
    cerr << "current step \t" << current_step << endl << endl;
    
    cerr << "total eqclasses \t" << eqclass_size.size() - 1 << endl;
  
    cerr << "eqc \t step \t above \t below" << endl;
    for(size_t i = 1; i < eqclass_size.size(); i++) {
      if(eqclass_size[i] != 0) {
        cerr << "e " << i << ":\t";
        cerr <<  eqclass_last_updated[i] << "\t";
        print_eqclass(eqclass_buddy_above[i]);
        cerr << "\t";
        print_eqclass(eqclass_buddy_below[i]);
        cerr << "\t";
        if(eqclass_buddy_above[i] == SINK_CLASS && eqclass_buddy_below[i] == SINK_CLASS) {
          cerr << "singleton";
        } else {
          cerr << "nonsingleton";
        }
        if(i == newest_eqclass) { cerr << "\t newest eqclass"; }
        cerr << endl;
      }
    }
    if(eqclass_size[newest_eqclass] == 0) {
      cerr << "e " << newest_eqclass << ":\t" << eqclass_last_updated[newest_eqclass] << "\t";
      print_eqclass(eqclass_buddy_above[newest_eqclass]);
      cerr << "\t";
      print_eqclass(eqclass_buddy_below[newest_eqclass]);
      cerr << "\t";
      if(eqclass_buddy_above[newest_eqclass] == SINK_CLASS && eqclass_buddy_below[newest_eqclass] == SINK_CLASS) {
        cerr << "singleton";
      } else {
        cerr << "nonsingleton";
      }
      cerr << "\t newest eqclass, empty";
      cerr << endl;
    }
    cerr << endl << "empty eqclass indices \t";
    for(size_t i = 0; i < empty_eqclass_indices.size(); i++) {
      cerr << empty_eqclass_indices[i] << "\t";
    }
    cerr << endl << endl;
    map_history.dump_state(current_step);
  }


  void mapHistory::dump_state(step_t marker) const {
    if(marker == CLEARED) {
      cerr << "argument was CLEARED" << endl;
    } else if(marker == PAST_FIRST) {
      cerr << "argument was PAST_FIRST" << endl;
    } else if(marker == NO_REP) {
      cerr << "argument was NO_REP" << endl;
    }
    
    // cerr << "n_eqclasses" << endl;
    // for(size_t i = 0; i < n_eqclasses.size(); i++) {
    //   if(i == marker) { cerr << "["; }
    //   cerr << n_eqclasses[i];
    //   if(i == marker) { cerr << "]"; }
    //   cerr << "\t";
    // }
    // cerr << endl;
    cerr << "step: \t";
    for(size_t i = 0; i < rep_eqclasses.size(); i++) {
      cerr << i << "\t";
    }
    cerr << "..." << endl;

    cerr << "reps:" << "\t";
    for(size_t i = 0; i < rep_eqclasses.size(); i++) {
      if(i == marker) { cerr << "["; }
      cerr << "e ";
      if(rep_eqclasses[i] == NO_REP) {
        cerr << "NR";
      } else {
        cerr << rep_eqclasses[i];
      }
      if(i == marker) { cerr << "]"; }
      cerr << "\t";
    }
    cerr << endl;
    cerr << "prev:" << "\t";
    for(size_t i = 0; i < previous.size(); i++) {
      if(i == marker) { cerr << "["; }
      if(previous[i] == CLEARED) {
        cerr << "CLRD";
      } else if(previous[i] == PAST_FIRST) {
        cerr << "PST1";
      } else {
        cerr << previous[i];
      }
      if(i == marker) { cerr << "]"; }
      cerr << "\t";
    }
    cerr << endl;
    cerr << "next:" << "\t";
    for(size_t i = 0; i < next.size(); i++) {
      if(i == marker) { cerr << "["; }
      if(next[i] == CLEARED) {
        cerr << "CLRD";
      } else {
        cerr << next[i];
      }
      if(i == marker) { cerr << "]"; }
      cerr << "\t";
    }
    cerr << endl;
  }
  
  bool mapHistory::validate_step(step_t i) const {
    if(i == CLEARED) {
      cerr << "invalid function argument: CLEARED" << endl;
    } else if(i == PAST_FIRST) {
      cerr << "invalid function argument: PAST_FIRST" << endl;
    } else if(i == NO_REP) {
      cerr << "invalid function argument: NO_REP" << endl;
    } else if(i >= size()) {
      cerr << "invalid function argument: out of bounds" << endl;
    }
    if(i >= size()) {
      dump_state(i);
      return false;
    }
    return true;
  }

  bool delayedEvalMap::ensure_unique(eqclass_t eqclass, step_t allowed_step) const {
    for(size_t i = 0; i <= newest_eqclass; i++) {
      if(!ensure_deleted(eqclass) && (eqclass_last_updated[i] != allowed_step)) {
        return false;
      }
    }
    for(size_t i = 0; i <= current_step; i++) {
      if(i != allowed_step && map_history.rep_eqclass(i) == eqclass) {
        return false;
      }
    }
    return true;
  }
  
  bool delayedEvalMap::ensure_deleted(eqclass_t eqclass) const {
    for(size_t i = 0; i < empty_eqclass_indices.size(); i++) {
      if(empty_eqclass_indices[i] == eqclass) {
        return true;
      }
    }
    return false;
  }
  
  bool delayedEvalMap::validate_eqclass(eqclass_t eqclass) const {
    if(eqclass == SINK_CLASS) {
      cerr << "invalid eqclass : SINK_CLASS" << endl;
      return false;
    }
    if(eqclass >= eqclass_size.size()) {
      cerr << "invalid eqclass : out of bounds" << endl;
      return false;
    }
    for(size_t j = 0; j < empty_eqclass_indices.size(); j++) {
      if(eqclass == empty_eqclass_indices[j]) {
        cerr << "invalid eqclass : deleted " << endl;
        return false;
      }
    }
    return true;
  }
  
  size_t delayedEvalMap::count_eqclasses_at_step(step_t i) const {
    eqclass_t eqclass = map_history.rep_eqclass(i);
    assert(validate_eqclass(eqclass));
    assert(eqclass_buddy_above[eqclass] == SINK_CLASS);
    size_t count = 1;
    while(eqclass_buddy_below[eqclass] != SINK_CLASS) {
      eqclass = eqclass_buddy_below[eqclass];
      count++;
    }
    return count;
  }
  
  bool delayedEvalMap::check_steps_valid() const {
    return check_steps_valid(vector<step_t>(0));
  }
  
  bool delayedEvalMap::check_steps_valid(vector<step_t> steps) const {
    vector<size_t> eqclass_counts_by_step(map_history.size(), 0);
    vector<size_t> rep_eqclass_counts(eqclass_size.size(), 0);
    vector<step_t> active_steps;
    step_t rightmost_first = 0;
    for(size_t i = map_history.size() - 1; i > 0; --i) {
      if(map_history.previous_step(i) == mapHistory::PAST_FIRST) {
        rightmost_first = i;
        break;
      }
    }
    assert(rightmost_first == map_history.leftmost_step());
    
    size_t total_count = 0;
    size_t steps_taken = 0;
    step_t it = current_step;
    while(it != mapHistory::PAST_FIRST) {
      active_steps.push_back(it);
      eqclass_t eqc_it = map_history.rep_eqclass(it);
      assert(eqclass_last_updated[map_history.rep_eqclass(it)] == it);
      assert(validate_eqclass(eqc_it));
      ++rep_eqclass_counts[eqc_it];
      eqclass_counts_by_step[it] = count_eqclasses_at_step(it);
      total_count += eqclass_counts_by_step[it];
      assert(it > map_history.previous_step(it) || map_history.previous_step(it) == mapHistory::PAST_FIRST);
      assert(it >= rightmost_first);
      it = map_history.previous_step(it);
      ++steps_taken;
      assert(steps_taken <= map_history.size());
    }
    assert(total_count == eqclass_size.size() - 1 - empty_eqclass_indices.size());
    for(size_t i = 0; i < rep_eqclass_counts.size(); i++) {
      assert(rep_eqclass_counts[i] < 2);
    }
    for(size_t i = 0; i < steps.size(); i++) {
      if(eqclass_counts_by_step[steps[i]] == 0) {
        cerr << "trying to update empty step" << endl;
        map_history.dump_state(steps[i]);
        return false;
      }
    }
    return true;
  }
#endif

bool delayedEvalMap::step_has_singleton_eqclass(step_t i) const {
#ifdef DEBUG
  assert(map_history.validate_step(i));
  assert(!((eqclass_buddy_below[map_history.rep_eqclass(i)] == SINK_CLASS) && (eqclass_buddy_above[map_history.rep_eqclass(i)] != SINK_CLASS)));
#endif
  return eqclass_buddy_below[map_history.rep_eqclass(i)] == SINK_CLASS; 
}

bool delayedEvalMap::eqclass_is_alone_at_step(eqclass_t eqclass) const {
#ifdef DEBUG
  assert(validate_eqclass(eqclass));
#endif  
  return eqclass_buddy_below[eqclass] == SINK_CLASS && eqclass_buddy_above[eqclass] == SINK_CLASS;  
}

void delayedEvalMap::push_back_at_nonempty_step(eqclass_t eqclass, step_t i) {
#ifdef DEBUG
  assert(map_history.validate_step(i));
  assert(validate_eqclass(eqclass));
  assert(eqclass_last_updated[eqclass] != i);
#endif
  eqclass_t old_rep = map_history.rep_eqclass(i);
  eqclass_buddy_above[eqclass] = SINK_CLASS;
  eqclass_buddy_below[eqclass] = old_rep;
  eqclass_buddy_above[old_rep] = eqclass;
  map_history.set_rep_eqclass(i, eqclass);
  eqclass_last_updated[eqclass] = i;
}

void delayedEvalMap::delete_non_singleton_eqclass_from_step(eqclass_t eqclass, step_t i) {
#ifdef DEBUG
  assert(map_history.validate_step(i));
  assert(validate_eqclass(eqclass));
  assert(!eqclass_is_alone_at_step(eqclass));
#endif
  if(eqclass == map_history.rep_eqclass(i)) {
    map_history.set_rep_eqclass(i, eqclass_buddy_below[eqclass]);
  }
  eqclass_buddy_above[eqclass_buddy_below[eqclass]] = eqclass_buddy_above[eqclass];
  eqclass_buddy_below[eqclass_buddy_above[eqclass]] = eqclass_buddy_below[eqclass];
}

void delayedEvalMap::move_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i) {
#ifdef DEBUG
  assert(map_history.validate_step(new_i));
  assert(old_i < current_step);
  assert(validate_eqclass(eqclass));
  assert(eqclass_last_updated[eqclass] == old_i);
  assert(old_i != map_history.leftmost_step());
#endif  
  if(eqclass_is_alone_at_step(eqclass)) {
    map_history.omit_step(eqclass_last_updated[eqclass]);
    push_back_at_nonempty_step(eqclass, new_i);
  } else {
    move_non_singleton_eqclass_to_nonempty_step(eqclass, old_i, new_i);
  }
}

void delayedEvalMap::move_leftmost_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i) {
#ifdef DEBUG
  assert(map_history.validate_step(new_i));
  assert(old_i < current_step);
  assert(validate_eqclass(eqclass));
  assert(eqclass_last_updated[eqclass] == old_i);
  assert(map_history.previous_step(eqclass_last_updated[eqclass]) == mapHistory::PAST_FIRST);
#endif  
  if(eqclass_is_alone_at_step(eqclass)) {
    map_history.set_new_leftmost(map_history.next_step(eqclass_last_updated[eqclass]));
#ifdef DEBUG
    map_history.clear_rep_eqclass(old_i);
#endif
    push_back_at_nonempty_step(eqclass, new_i);
  } else {
    move_non_singleton_eqclass_to_nonempty_step(eqclass, old_i, new_i);
  }
}

void delayedEvalMap::move_non_singleton_eqclass_to_nonempty_step(eqclass_t eqclass, step_t old_i, step_t new_i) {
#ifdef DEBUG
  assert(map_history.validate_step(old_i));
  assert(map_history.validate_step(new_i));
  assert(validate_eqclass(eqclass));
  assert(eqclass_last_updated[eqclass] == old_i);
  assert(!eqclass_is_alone_at_step(eqclass));
#endif
  if(eqclass == map_history.rep_eqclass(old_i)) {
    map_history.set_rep_eqclass(old_i, eqclass_buddy_below[eqclass]);
  }
  eqclass_buddy_above[eqclass_buddy_below[eqclass]] = eqclass_buddy_above[eqclass];
  eqclass_buddy_below[eqclass_buddy_above[eqclass]] = eqclass_buddy_below[eqclass];
  push_back_at_nonempty_step(eqclass, new_i);
}

void delayedEvalMap::add_empty_eqclass_to_empty_step(const DPUpdateMap& map) {
  if(empty_eqclass_indices.size() == 0) {
    newest_eqclass = eqclass_to_map.size();
    eqclass_to_map.push_back(map);
    eqclass_size.push_back(0);
    eqclass_last_updated.push_back(current_step);
    eqclass_buddy_above.push_back(SINK_CLASS);
    eqclass_buddy_below.push_back(SINK_CLASS);
  } else {
    newest_eqclass = empty_eqclass_indices.back();
#ifdef DEBUG
    assert(ensure_deleted(newest_eqclass));
    assert(ensure_unique(newest_eqclass, current_step));
#endif
    empty_eqclass_indices.pop_back();
    eqclass_to_map[newest_eqclass] = map;
    eqclass_size[newest_eqclass] = 0;
    eqclass_last_updated[newest_eqclass] = current_step;
    eqclass_buddy_above[newest_eqclass] = SINK_CLASS;
    eqclass_buddy_below[newest_eqclass] = SINK_CLASS;
  }
  map_history.assign_rep_eqclass_to_new_step(newest_eqclass);
  map_history.assign_map_to_new_step(map);
}

void delayedEvalMap::move_row_to_newest_eqclass(row_t row) {
#ifdef DEBUG
  assert(eqclass_last_updated[row_to_eqclass[row]] == current_step);
#endif
  decrement_eqclass(row_to_eqclass[row]);
  row_to_eqclass[row] = newest_eqclass;
  ++eqclass_size[newest_eqclass];
}

void delayedEvalMap::decrement_eqclass(eqclass_t eqclass) {
  if(eqclass_size[eqclass] == 1) {
    delete_eqclass_from_step(eqclass);
  } else {
    --eqclass_size[eqclass];
  }
}

void delayedEvalMap::delete_eqclass_from_step(eqclass_t eqclass) {
  size_t step = eqclass_last_updated[eqclass];
  if(map_history.rep_eqclass(step) == eqclass) {
    map_history.set_rep_eqclass(step, eqclass_buddy_below[eqclass]);
  }
  
  empty_eqclass_indices.push_back(eqclass);
  eqclass_buddy_above[eqclass_buddy_below[eqclass]] = eqclass_buddy_above[eqclass];
  eqclass_buddy_below[eqclass_buddy_above[eqclass]] = eqclass_buddy_below[eqclass];
#ifdef DEBUG
  if(map_history.rep_eqclass(step) == SINK_CLASS) {
    map_history.clear_rep_eqclass(step);
  }
  eqclass_size[eqclass] = 0;
#endif
}

vector<eqclass_t> delayedEvalMap::rows_to_eqclasses(const rowSet& rows) const {
  #ifdef DEBUG
  assert(!rows.empty());
  #endif
  
  vector<eqclass_t> to_return(0);
  if(rows.empty()) {
    return to_return;
  }
  vector<char> seen = vector<char>(eqclass_to_map.size(), 0);
  rowSet::const_iterator it = rows.begin();
  rowSet::const_iterator rows_end = rows.end();
  for(it; it != rows_end; ++it) {
    if(seen[row_to_eqclass[*it]] == 0) {
      to_return.push_back(row_to_eqclass[*it]);
    }
    ++seen[row_to_eqclass[*it]];
  }
  
#ifdef DEBUG
  assert(to_return.size() > 0);
  for(size_t i = 0; i < to_return.size(); i++) {
    assert(validate_eqclass(to_return[i]));
  }
#endif
  
  return to_return;
}

void delayedEvalMap::extend_by_new_step(const DPUpdateMap& map) {
  map_history.increment_step();
  ++current_step;
  add_empty_eqclass_to_empty_step(map);
#ifdef DEBUG
  assert(map_history.rep_eqclass(current_step) == newest_eqclass);  assert(step_has_singleton_eqclass(current_step));
  assert(eqclass_size[newest_eqclass] == 0);
  assert(validate_eqclass(newest_eqclass));
#endif
}

void delayedEvalMap::extend_value_only(const DPUpdateMap& map) {
  map_history.postcompose_rightmost_map(map);
}

void delayedEvalMap::update_evaluate_and_move_rows(const rowSet& active_rows, vector<double>& values, double minority_correction) {
  vector<eqclass_t> active_eqclasses = rows_to_eqclasses(active_rows);
  update_eqclasses(active_eqclasses);
  rowSet::const_iterator it = active_rows.begin();
  rowSet::const_iterator active_rows_end = active_rows.end();
  for(it; it != active_rows_end; ++it) {
    values[*it] = eqclass_to_map[row_to_eqclass[*it]].of(values[*it]) + minority_correction;
    move_row_to_newest_eqclass(*it);
  }
}

void delayedEvalMap::update_eqclasses(const vector<eqclass_t>& eqclasses) {
#ifdef DEBUG
  assert(check_steps_valid());
#endif
  step_t oldest_step = current_step;
  for(eqclass_t i = 0; i < eqclasses.size(); i++) {
#ifdef DEBUG
    assert(eqclass_last_updated[eqclasses[i]] < current_step);
    assert(map_history.validate_step(eqclass_last_updated[eqclasses[i]]));
#endif
    if(eqclass_last_updated[eqclasses[i]] < oldest_step) {
      oldest_step = eqclass_last_updated[eqclasses[i]];
    }
  }
  
  if(current_step == oldest_step) {
    return;
  }
  
#ifdef DEBUG
  assert(oldest_step >= map_history.leftmost_step());
  assert(check_steps_valid());
#endif
  
  step_t last_i = current_step;
  map_history.suffix(last_i) = map_history[current_step];
  
  for(step_t i = map_history.previous_step(last_i); last_i > oldest_step; ) {
    map_history.suffix(i) = map_history.suffix(last_i).compose(map_history[i]);
    last_i = i;
    i = map_history.previous_step(last_i);
  }
  
  last_i = current_step;
    
  for(step_t i = map_history.previous_step(current_step); i > oldest_step; ) {
    if(step_has_singleton_eqclass(i)) {
      eqclass_t i_rep_eqclass = map_history.rep_eqclass(i);
#ifdef DEBUG
      assert(validate_eqclass(i_rep_eqclass));
#endif
      // we skip clearing eqclass_buddy_above/below because this gets overwritten by the following call
      push_back_at_nonempty_step(i_rep_eqclass, current_step);
      eqclass_to_map[i_rep_eqclass] = map_history.suffix(last_i).of(eqclass_to_map[i_rep_eqclass]);
      map_history.omit_step(i);
    } else {
      last_i = i;
    }
    i = map_history.previous_step(last_i);
  }

#ifdef DEBUG
  assert(check_steps_valid());
#endif
  
  for(size_t j = 0; j < eqclasses.size(); j++) {
    eqclass_t eqclass = eqclasses[j];
    step_t last_update = eqclass_last_updated[eqclass];
    if(last_update != current_step) {
      eqclass_to_map.at(eqclass) = map_history.suffix(last_update).of(eqclass_to_map.at(eqclass));
      if(last_update != map_history.leftmost_step()) {
        move_eqclass_to_nonempty_step(eqclass, last_update, current_step);
      } else {
        move_leftmost_eqclass_to_nonempty_step(eqclass, last_update, current_step);
      }
    }
  }
  
#ifdef DEBUG
  assert(check_steps_valid());
#endif
}

// ---------------------------------------------------------------------------------------------------------------------

delayedEvalMap::delayedEvalMap() {
  
}

delayedEvalMap::delayedEvalMap(size_t rows) : 
	row_to_eqclass(vector<size_t>(rows, 1)), 
	eqclass_last_updated(vector<size_t>(2, 0)),
	newest_eqclass(1),
	eqclass_to_map(vector<DPUpdateMap>(2, DPUpdateMap::IDENTITY)),
	map_history(mapHistory(DPUpdateMap::IDENTITY)),
  eqclass_buddy_above(2, SINK_CLASS),
  eqclass_buddy_below(2, SINK_CLASS) {
  eqclass_size = {0, rows};
}

void delayedEvalMap::move_all_rows_to_newest_eqclass() {
  for(row_t row = 0; row < row_to_eqclass.size(); row++) {
    move_row_to_newest_eqclass(row);
  }
}

void delayedEvalMap::hard_update_all() {
  vector<size_t> non_empty_eqclasses;
  for(size_t i = 1; i < eqclass_size.size(); i++) {
    if(eqclass_size[i] != 0) {
      non_empty_eqclasses.push_back(i);
    }
  }
  update_eqclasses(non_empty_eqclasses);
}

const DPUpdateMap& delayedEvalMap::get_map(row_t row) const {
  return eqclass_to_map[row_to_eqclass[row]];
}

double delayedEvalMap::evaluate(size_t row, double value) const {
  return eqclass_to_map[row_to_eqclass[row]].of(value);
}