#include "row_set.hpp"
#include "delay_multiplier.hpp"
#include "math.hpp"
#include <iostream>

using namespace std;

// TODO TODO TODO
// deleting last singleton doesn't currently account for composition of maps correctly! Do we fuse with below?
// TODO TODO  TODO




// ---------------------------------------------------------------------------------------------------------------------
// -- mapHistory -------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------

const step_t mapHistory::CLEARED = SIZE_MAX - 4;
const step_t mapHistory::PAST_FIRST = SIZE_MAX - 3;
const eqclass_t mapHistory::NO_REP = SIZE_MAX - 2;

const eqclass_t delayedEvalMap::INACTIVE_EQCLASS = SIZE_MAX - 1;
const eqclass_t delayedEvalMap::NO_NEIGHBOUR = SIZE_MAX;

mapHistory::mapHistory() {
}

mapHistory::mapHistory(const DPUpdateMap& map) {
	elements = {map};
  previous = {PAST_FIRST, 0};
  next = {1, 2};
  suffixes = {DPUpdateMap::IDENTITY};
  rep_eqclasses = {0};
  n_eqclasses = {1};
}

mapHistory::mapHistory(const mapHistory& other) {
  elements = other.elements;
  previous = other.previous;
  next = other.next;
  suffixes = other.suffixes;
}

void mapHistory::reserve(size_t length) {
  elements.reserve(length);
  suffixes.reserve(length);
  previous.reserve(length);
  next.reserve(length);
  rep_eqclasses.reserve(length);
}

step_t mapHistory::max_index() const {
  return elements.size() - 1;
}

void mapHistory::push_back(const DPUpdateMap& map) {
	elements.push_back(map);
  previous.push_back(previous.size() - 1);
  next.push_back(next.size() + 1);
  suffixes.push_back(DPUpdateMap::IDENTITY);
  rep_eqclasses.push_back(NO_REP);
  n_eqclasses.push_back(0);
}

const DPUpdateMap& mapHistory::operator[](step_t i) const {
  #ifdef DEBUG
    return elements.at(i);
  #endif
	return elements[i];
}

DPUpdateMap& mapHistory::operator[](step_t i) {
  #ifdef DEBUG
    return elements.at(i);
  #endif
	return elements[i];
}

const DPUpdateMap& mapHistory::back() const {
	return elements.back();
}

DPUpdateMap& mapHistory::back() {
	return elements.back();
}

const DPUpdateMap& mapHistory::suffix(step_t i) const {
  #ifdef DEBUG
    const DPUpdateMap& this_suffix = suffixes.at(i);
    if(previous.at(i) == CLEARED) {
      throw erased_error("Wanted suffix of cleared step "+i);
    }
    return this_suffix;
  #endif
  return suffixes[i];
}

DPUpdateMap& mapHistory::suffix(step_t i) {
  #ifdef DEBUG
    DPUpdateMap& this_suffix = suffixes.at(i);
    if(previous.at(i) == CLEARED) {
      throw erased_error("Wanted suffix of cleared step "+i);
    }
  #endif
  return suffixes[i];
}

step_t mapHistory::previous_step(step_t i) const {
  #ifdef DEBUG
    step_t this_previous = previous.at(i);
    if(this_previous == CLEARED) {
      throw erased_error("Wanted previous of cleared step "+i);
    }
    if(this_previous == PAST_FIRST) {
      throw erased_error("Wanted previous of first step "+i);
    }
  #endif
  return previous[i];
}

step_t mapHistory::next_step(step_t i) const {
  #ifdef DEBUG
    step_t this_next = next.at(i);
    if(this_next == CLEARED) {
      throw erased_error("Wanted next of cleared step "+i);
    }
    if(this_next >= next.size() - 1) {
      throw runtime_error("Wanted out of range next step");
    }
  #endif
  return next[i];
}

bool mapHistory::step_cleared(step_t i) const {
  #ifdef DEBUG
    return next.at(i) == CLEARED;
  #else
    return next[i] == CLEARED;
  #endif
}

void mapHistory::fuse_prev(size_t i) {
  #ifdef DEBUG
    step_t old_previous = previous.at(i);
    step_t old_next = next.at(i);
    if(old_previous == CLEARED) {
      throw erased_error("Tried to fuse previous of cleared step "+i);
    }
    if(old_previous == PAST_FIRST) {
      throw erased_error("Wanted previous of first step "+i);
    }
    suffixes.at(i) = suffixes.at(old_previous);
    previous.at(i) = previous_step(old_previous);
    next.at(previous.at(i)) = i;
    elements.at(i) = elements.at(i).of(elements.at(old_previous));
  #else
    step_t old_previous = previous[i];
    suffixes[i] = suffixes[old_previous];
    previous[i] = previous[old_previous];
    next[previous[i]] = i;
    elements[i] = elements[i].of(elements[old_previous]);
  #endif
  previous[old_previous] = CLEARED;
  next[old_previous] = CLEARED;
  rep_eqclasses[old_previous] = NO_REP;
  n_eqclasses[old_previous] = 0;
}

size_t mapHistory::size() const {
	return elements.size();
}

const vector<DPUpdateMap>& mapHistory::get_elements() const {
  return elements;
}

void mapHistory::set_rep_eqclass(step_t i, eqclass_t eqclass) {
  #ifdef DEBUG
    for(size_t j = 0; j < rep_eqclasses.size(); j++) {
      if(rep_eqclasses[j] == eqclass) {
        throw runtime_error("set eqclass as representative of two sites");
      }
    }
    rep_eqclasses.at(i) = eqclass;
  #else
    rep_eqclasses[i] = eqclass;
  #endif
}

eqclass_t mapHistory::rep_eqclass(step_t i) const {
  #ifdef DEBUG
    return rep_eqclasses.at(i);
  #else
    return rep_eqclasses[i];
  #endif
}

void mapHistory::clear_end_singleton() {
  step_t i = max_index();
  step_t past_end = i + 1;
  step_t this_previous = previous.at(i);
  #ifdef DEBUG
    rep_eqclasses.at(i) = NO_REP;
    next.at(i) = CLEARED;
    previous.at(i) = CLEARED;
    n_eqclasses.at(i) = 0;
    if(this_previous != PAST_FIRST && this_previous != NO_REP && this_previous != CLEARED) {
      next.at(this_previous) = past_end;
    }
    previous.at(past_end) = this_previous;
  #else
    rep_eqclasses[i] = NO_REP;
    next[i] = CLEARED;
    previous[i] = CLEARED;
    n_eqclasses[i] = 0;
    if(this_previous != PAST_FIRST && this_previous != NO_REP && this_previous != CLEARED) {
      next[this_previous] = past_end;
    }
    previous[past_end] = this_previous;
  #endif
}

void mapHistory::clear_singleton(step_t i) {
  #ifdef DEBUG
    if(n_eqclasses.at(i) != 1) {
      throw runtime_error("called clear_singleton on non-singleton");
    }
    if(previous.at(i) != PAST_FIRST) {
      fuse_prev(next.at(i));
    } else {
      step_t this_next = next.at(i);
      previous.at(this_next) = PAST_FIRST;
      n_eqclasses.at(i) = 0;
    }
    rep_eqclasses.at(i) = NO_REP;
    next.at(i) = CLEARED;
    previous.at(i) = CLEARED;
  #else
    rep_eqclasses[i] = NO_REP;
    next[i] = CLEARED;
    prev[i] = CLEARED;
  #endif
}

inline bool delayedEvalMap::is_singleton(eqclass_t eqclass) const {
  #ifdef DEBUG
    if(eqclass == NO_NEIGHBOUR) {
      throw runtime_error("called is_singleton on NO_NEIGHBOUR");
    } else if(eqclass == mapHistory::NO_REP) {
      throw runtime_error("called is_singleton onmapHistory::NO_REP");
    } else if(eqclass == INACTIVE_EQCLASS) {
      throw runtime_error("called is_singleton on INACTIVE_EQCLASS");
    }
    return map_history.n_eqclasses[eqclass_last_updated[eqclass]] == 1;
    // return (eqclass_buddy_above.at(eqclass) == NO_NEIGHBOUR) && (eqclass_buddy_below.at(eqclass) == NO_NEIGHBOUR);
  #endif
  return (eqclass_buddy_above[eqclass] == NO_NEIGHBOUR) && (eqclass_buddy_below[eqclass] == NO_NEIGHBOUR);
}

inline bool delayedEvalMap::is_front(eqclass_t eqclass) const {
  #ifdef DEBUG
    if(eqclass == NO_NEIGHBOUR) {
      throw runtime_error("called is_front on NO_NEIGHBOUR");
    } else if(eqclass == mapHistory::NO_REP) {
      throw runtime_error("called is_front onmapHistory::NO_REP");
    } else if(eqclass == INACTIVE_EQCLASS) {
      throw runtime_error("called is_front on INACTIVE_EQCLASS");
    }
  #endif
  return (eqclass_buddy_above[eqclass] == NO_NEIGHBOUR);
}

void delayedEvalMap::push_back_at_current(eqclass_t eqclass) {
  #ifdef DEBUG
    if(eqclass == NO_NEIGHBOUR) {
      throw runtime_error("called push_back_at_current on NO_NEIGHBOUR");
    } else if(eqclass == mapHistory::NO_REP) {
      throw runtime_error("called push_back_at_current onmapHistory::NO_REP");
    } else if(eqclass == INACTIVE_EQCLASS) {
      throw runtime_error("called push_back_at_current on INACTIVE_EQCLASS");
    }
    if(map_history.rep_eqclass(current_site) == mapHistory::NO_REP) {
      eqclass_buddy_above.at(eqclass) = NO_NEIGHBOUR;
      eqclass_buddy_below.at(eqclass) = NO_NEIGHBOUR;
      map_history.set_rep_eqclass(current_site, eqclass);
    } else {
      eqclass_t current_rep = map_history.rep_eqclass(current_site);
      eqclass_buddy_below.at(eqclass) = current_rep;
      eqclass_buddy_above.at(current_rep) = eqclass;
      eqclass_buddy_above.at(eqclass) = NO_NEIGHBOUR;
      map_history.set_rep_eqclass(current_site, eqclass);
    }
  #else
    if(map_history.rep_eqclass(current_site) == mapHistory::NO_REP) {
      eqclass_buddy_above[eqclass] = NO_NEIGHBOUR;
      eqclass_buddy_below[eqclass] = NO_NEIGHBOUR;
      map_history.set_rep_eqclass(current_site, eqclass);
    } else {
      eqclass_t current_rep = map_history.rep_eqclass(current_site);
      eqclass_buddy_below[eqclass] = current_rep;
      eqclass_buddy_above[current_rep] = eqclass;
      eqclass_buddy_above[eqclass] = NO_NEIGHBOUR;
      map_history.set_rep_eqclass(current_site, eqclass);
    }
  #endif
  ++map_history.n_eqclasses[current_site];
}

void delayedEvalMap::delete_at_site(eqclass_t eqclass) {
  size_t departing_site = eqclass_last_updated.at(eqclass);
  #ifdef DEBUG
    if(eqclass == NO_NEIGHBOUR) {
      throw runtime_error("called delete_at_site on NO_NEIGHBOUR");
    } else if(eqclass == mapHistory::NO_REP) {
      throw runtime_error("called delete_at_site onmapHistory::NO_REP");
    } else if(eqclass == INACTIVE_EQCLASS) {
      throw runtime_error("called delete_at_site on INACTIVE_EQCLASS");
    }
    if(is_singleton(eqclass)) {
      step_t this_step = eqclass_last_updated.at(eqclass); 
      eqclass_buddy_above.at(eqclass) == INACTIVE_EQCLASS;
      eqclass_buddy_below.at(eqclass) == INACTIVE_EQCLASS;
      if(this_step == map_history.max_index()) {
        map_history.clear_end_singleton();
      } else {
        map_history.clear_singleton(eqclass_last_updated.at(eqclass));
      }
    } else {
      eqclass_t old_above = eqclass_buddy_above.at(eqclass);
      eqclass_t old_below = eqclass_buddy_below.at(eqclass);
      if(old_above != NO_NEIGHBOUR) {
        eqclass_buddy_below.at(old_above) = old_below;
      } else {
        // safe to not check because the case that both above and below are NO_NEIGHBOUR is caught by is_singleton(eqclass)
        map_history.set_rep_eqclass(eqclass_last_updated.at(eqclass), old_below);
      }
      if(old_below != NO_NEIGHBOUR) {
        eqclass_buddy_above.at(old_below) = old_above;
      }
      eqclass_buddy_above.at(eqclass) == INACTIVE_EQCLASS;
      eqclass_buddy_below.at(eqclass) == INACTIVE_EQCLASS;
      // if(map_history.n_eqclasses[departing_site] == 0) {
      //   throw runtime_error("double cleared site");
      // }
      --map_history.n_eqclasses[departing_site];
    }
  #else
    if(is_singleton(eqclass)) {
      step_t this_step = eqclass_last_updated[eqclass]; 
      if(this_step == map_history.max_index()) {
        map_history.clear_end_singleton();
      } else {
        map_history.clear_singleton(eqclass_last_updated.at(eqclass));
      }
    } else {
      eqclass_t old_above = eqclass_buddy_above.at(eqclass);
      eqclass_t old_below = eqclass_buddy_below.at(eqclass);
      if(old_above != NO_NEIGHBOUR) {
        eqclass_buddy_below.at(old_above) = old_below;
      } else {
        // safe to not check because the case that both above and below are NO_NEIGHBOUR is caught by is_singleton(eqclass)
        map_history.set_rep_eqclass(eqclass_last_updated.at(eqclass), old_below);
      }
      if(old_below != NO_NEIGHBOUR) {
        eqclass_buddy_above.at(old_below) = old_above;
      }
    }
  #endif
}

// ---------------------------------------------------------------------------------------------------------------------

delayedEvalMap::delayedEvalMap() {
  
}

delayedEvalMap::delayedEvalMap(size_t rows) : 
	row_to_eqclass(vector<size_t>(rows, 0)), 
	eqclass_size(vector<size_t>(1, rows)),
	eqclass_last_updated(vector<size_t>(1, 0)),
	newest_eqclass(0),
	eqclass_to_map(vector<DPUpdateMap>(1, DPUpdateMap::IDENTITY)),
	map_history(mapHistory(DPUpdateMap::IDENTITY)),
  eqclass_buddy_above(1, NO_NEIGHBOUR),
  eqclass_buddy_below(1, NO_NEIGHBOUR) {
}

void delayedEvalMap::add_identity_eqclass() {
  add_eqclass(DPUpdateMap::IDENTITY);
  return;
}

delayedEvalMap::delayedEvalMap(const delayedEvalMap &other) {
	current_site = other.current_site;
	row_to_eqclass = other.row_to_eqclass;
	eqclass_last_updated = other.eqclass_last_updated;
  size_t oldest_seen = current_site;
  for(size_t i = 0; i < eqclass_last_updated.size(); i++) {
    if(eqclass_last_updated[i] < oldest_seen) {
      oldest_seen = eqclass_last_updated[i];
    }
  }
  map_history = mapHistory(other.map_history);
	eqclass_to_map = other.eqclass_to_map;
	eqclass_size = other.eqclass_size;
	empty_eqclass_indices = other.empty_eqclass_indices;
}

void delayedEvalMap::assign_row_to_newest_eqclass(size_t row) {
  row_to_eqclass[row] = newest_eqclass;
  eqclass_size[newest_eqclass]++;
  return;
}

void delayedEvalMap::hard_clear_all() {
  for(int i = 0; i < eqclass_to_map.size(); i++) {
    delete_eqclass(i);
  }
  add_identity_eqclass();
  for(int i = 0; i < row_to_eqclass.size(); i++) {
    assign_row_to_newest_eqclass(i);
  }
  return;
}

void delayedEvalMap::hard_update_all() {
  vector<size_t> non_empty_eqclasses;
  
  for(int i = 0; i < eqclass_size.size(); i++) {
    if(eqclass_size[i] != 0) {
      non_empty_eqclasses.push_back(i);
    }
  }
  update_maps(non_empty_eqclasses);
  return;
}

vector<size_t> delayedEvalMap::rows_to_eqclasses(const rowSet& rows) const {
  vector<size_t> to_return(0);
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
    seen[row_to_eqclass[*it]] = 1;
  }
  return to_return;
}

// implement bounds checking on all map_history

void delayedEvalMap::update_maps(const vector<size_t>& eqclasses) {
  step_t least_up_to_date = current_site;
  for(eqclass_t i = 0; i < eqclasses.size(); i++) {
    if(eqclass_last_updated[eqclasses[i]] < least_up_to_date) {
      least_up_to_date = eqclass_last_updated[eqclasses[i]];
    }
  }
  if(current_site != least_up_to_date) {
    // ONLY DIFFERENCES AFTER THIS POINT
    step_t last_i = current_site;
    
    map_history.suffix(current_site) = map_history[current_site];
    
    for(step_t i = current_site; i > least_up_to_date; ) {
      i = map_history.previous_step(i);
      map_history.suffix(i) = map_history.suffix(last_i).compose(map_history[i]);
      last_i = i;
    }
    
    last_i = current_site;
    
    for(step_t i = current_site - 1; i > least_up_to_date; ) {
      eqclass_t eqclass = map_history.rep_eqclass(i); 
      if(eqclass == mapHistory::NO_REP) {
        map_history.fuse_prev(last_i);
        i = map_history.previous_step(last_i);
      } else if(is_singleton(eqclass)) {
        // HAPPENING IN HERE
        eqclass_to_map[eqclass] = map_history.suffix(last_i).of(eqclass_to_map[eqclass]);
        delete_at_site(eqclass);
        push_back_at_current(eqclass);
        eqclass_last_updated[eqclass] = current_site;
        map_history.fuse_prev(last_i);
        i = map_history.previous_step(last_i);
      } else {
        last_i = i;
        i = map_history.previous_step(last_i);
      }
    }    
    
    for(size_t i = 0; i < eqclasses.size(); i++) {
      #ifdef DEBUG
      eqclass_t eqclass = eqclasses.at(i);
      if(eqclass_last_updated.at(eqclass) != current_site && !map_history.step_cleared(eqclass_last_updated.at(eqclass))) {
        // j is the eqclass's index in the suffix-vector
        size_t j = map_history.next_step(eqclass_last_updated[eqclass]);
        eqclass_to_map.at(eqclass) = map_history.suffix(j).of(eqclass_to_map.at(eqclass));
        delete_at_site(eqclass);
        push_back_at_current(eqclass);
        eqclass_last_updated.at(eqclass) = current_site;
      }
      #else
      eqclass = eqclasses[i];
      if(eqclass_last_updated[eqclass] != current_site &&  map_history.next_step(eqclass_last_updated.at(eqclass)) != mapHistory::CLEARED) {
        // j is the eqclass's index in the suffix-vector
        size_t j = map_history.next_step(eqclass_last_updated[eqclass]);
        eqclass_to_map[eqclass] = map_history.suffix(j).of(eqclass_to_map[eqclass]);
        delete_at_site(eqclass);
        push_back_at_current(eqclass);
        eqclass_last_updated[eqclass] = current_site;
      }
      #endif
    }
  }
  return;
}

void delayedEvalMap::delete_eqclass(size_t eqclass) {
  delete_at_site(eqclass);
  eqclass_size[eqclass] = 0;
  eqclass_last_updated[eqclass] = current_site;
  empty_eqclass_indices.push_back(eqclass);
  return;
}

void delayedEvalMap::decrement_eqclass(size_t eqclass) {
  if(eqclass_size[eqclass] == 1) {
    delete_eqclass(eqclass);
  } else {
    --eqclass_size[eqclass];
  }
  return;
}

void delayedEvalMap::remove_row_from_eqclass(size_t row) {
  decrement_eqclass(row_to_eqclass[row]);
  // unassigned row is given max possible eqclass index + 1 to ensure that
  // accessing it will throw an error
  row_to_eqclass[row] = row_to_eqclass.size();
  return;
}

void delayedEvalMap::add_eqclass(const DPUpdateMap& map) {
  if(empty_eqclass_indices.size() == 0) {
    newest_eqclass = eqclass_to_map.size();
    eqclass_to_map.push_back(map);
    eqclass_size.push_back(0);
    eqclass_last_updated.push_back(current_site);
    eqclass_buddy_above.push_back(NO_NEIGHBOUR);
    eqclass_buddy_below.push_back(NO_NEIGHBOUR);
    push_back_at_current(newest_eqclass);
    return;
  } else {
    newest_eqclass = empty_eqclass_indices.back();
    empty_eqclass_indices.pop_back();
    eqclass_to_map[newest_eqclass] = map;
    eqclass_size[newest_eqclass] = 0;
    eqclass_last_updated[newest_eqclass] = current_site;
    push_back_at_current(newest_eqclass);
    return;
  }
}

double delayedEvalMap::get_constant(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]].constant;
}

double delayedEvalMap::get_coefficient(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]].coefficient;
}

const DPUpdateMap& delayedEvalMap::get_map(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]];
}

const vector<DPUpdateMap>& delayedEvalMap::get_maps() const {
  return eqclass_to_map;
}

vector<DPUpdateMap>& delayedEvalMap::get_maps() {
  return eqclass_to_map;
}

const vector<size_t>& delayedEvalMap::get_map_indices() const {
  return row_to_eqclass;
}

void delayedEvalMap::stage_map_for_span(const DPUpdateMap& span_map) {
  stage_map_for_site(span_map);
  return;
}

void delayedEvalMap::stage_map_for_site(const DPUpdateMap& site_map) {
  current_site++;
  map_history.push_back(site_map);
  return;
}

size_t delayedEvalMap::last_update(size_t row) const {
  if(row_to_eqclass[row] != row_to_eqclass.size()) {
    return eqclass_last_updated[row_to_eqclass[row]];
  } else {
    return current_site;
  }
}

const vector<DPUpdateMap>& delayedEvalMap::get_map_history() const {
  return map_history.get_elements();
}

void delayedEvalMap::reset_rows(const rowSet& rows) {
  rowSet::const_iterator it = rows.begin();
  rowSet::const_iterator rows_end = rows.end();
  for(it; it != rows_end; ++it) {
    remove_row_from_eqclass(*it);
  }
  add_identity_eqclass();
  
  it = rows.begin();
  for(it; it != rows_end; ++it) {
    assign_row_to_newest_eqclass(*it);
  }
}

void delayedEvalMap::update_active_rows(const rowSet& active_rows) {
  update_maps(rows_to_eqclasses(active_rows));
}

size_t delayedEvalMap::number_of_eqclasses() const {
  return eqclass_size.size() - empty_eqclass_indices.size();
}

size_t delayedEvalMap::row_updated_to(size_t row) const {
  return eqclass_last_updated[row_to_eqclass[row]];
}

size_t delayedEvalMap::get_current_site() const {
  return current_site;
}

size_t delayedEvalMap::get_eqclass(size_t row) const {
  return row_to_eqclass[row];
}

double delayedEvalMap::evaluate(size_t row, double value) const {
  return eqclass_to_map[row_to_eqclass[row]].of(value);
}
