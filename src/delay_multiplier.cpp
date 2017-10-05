#include "row_set.hpp"
#include "delay_multiplier.hpp"
#include "math.hpp"
#include <iostream>

using namespace std;

void mapHistory::push_back(const DPUpdateMap& map) {
	elements.push_back(map);
}

size_t mapHistory::size() const {
	return elements.size();
}

mapHistory::mapHistory() {
  start = 0;
}

mapHistory::mapHistory(const DPUpdateMap& map, size_t start) : start(start) {
	elements = {map};
}


mapHistory::mapHistory(const mapHistory& other) {
	start = other.start;
  elements = other.elements;
}

mapHistory::mapHistory(const mapHistory& other, size_t new_start) {
	start = new_start;
	size_t offset = new_start - other.start;
	elements = vector<DPUpdateMap>(other.elements.begin() + offset, other.elements.end());	
}

DPUpdateMap& mapHistory::operator[](size_t i) {
	return elements[i - start];
}

DPUpdateMap& mapHistory::back() {
	return elements.back();
}

const vector<DPUpdateMap>& mapHistory::get_elements() const {
  return elements;
}

delayMap::delayMap() {
  
}

delayMap::delayMap(size_t rows, size_t start) : 
            dM_start(start) {
  maps_by_site = mapHistory(DPUpdateMap(0), start);
  current_site = start;
  add_identity_map();
  // After calling the above; maps_by_slot is a singleton containing the
  // identity map. updated_to has been set for this slot; count has not
  counts = {rows};
  // All rows correspond to this slot
  slots_by_row = vector<size_t>(rows, 0);
}

void delayMap::add_identity_map() {
  add_map(DPUpdateMap(0));
  return;
}

delayMap::delayMap(const delayMap &other) {
	added_span = other.added_span;
	updated_maps = other.updated_maps;
	// dM_start = other.dM_start;
	current_site = other.current_site;
	slots_by_row = other.slots_by_row;
	updated_to = other.updated_to;
  size_t oldest_seen = current_site;
  for(size_t i = 0; i < updated_to.size(); i++) {
    if(updated_to[i] < oldest_seen) {
      oldest_seen = updated_to[i];
    }
  }
  dM_start = oldest_seen;
  maps_by_site = mapHistory(other.maps_by_site, oldest_seen);
	maps_by_slot = other.maps_by_slot;
	counts = other.counts;
	empty_map_slots = other.empty_map_slots;
}

void delayMap::assign_row_to_newest_index(size_t row) {
  //TODO: complain if slots_by_row[row] != |H|
  slots_by_row[row] = newest_index;
  counts[newest_index]++;
  return;
}

void delayMap::hard_clear_all() {
  for(int i = 0; i < maps_by_slot.size(); i++) {
    delete_slot(i);
  }
  add_identity_map();
  for(int i = 0; i < slots_by_row.size(); i++) {
    assign_row_to_newest_index(i);
  }
  return;
}

void delayMap::hard_update_all() {
  vector<size_t> non_empty_slots;
  
  for(int i = 0; i < counts.size(); i++) {
    if(counts[i] != 0) {
      non_empty_slots.push_back(i);
    }
  }
  update_maps(non_empty_slots);
  return;
}

vector<size_t> delayMap::rows_to_slots(const vector<size_t>& rows) const {
  vector<bool> seen = vector<bool>(maps_by_slot.size(), false);
  vector<size_t> to_return;
  for(int i = 0; i < rows.size(); i++) {
    if(!(seen[slots_by_row[rows[i]]])) {
      to_return.push_back(slots_by_row[rows[i]]);
    }
    seen[slots_by_row[rows[i]]] = true;
  }
  return to_return;
}

vector<size_t> delayMap::rows_to_slots(const rowSet& rows) const {
  vector<bool> seen = vector<bool>(maps_by_slot.size(), false);
  vector<size_t> to_return;
  for(int i = 0; i < rows.size(); i++) {
    if(!(seen[slots_by_row[rows[i]]])) {
      to_return.push_back(slots_by_row[rows[i]]);
    }
    seen[slots_by_row[rows[i]]] = true;
  }
  return to_return;
}

vector<bool> delayMap::rows_to_slotmask(const rowSet& rows) const {
  vector<bool> to_return = vector<bool>(maps_by_slot.size(), false);
  for(int i = 0; i < rows.size(); i++) {
    to_return[slots_by_row[rows[i]]] = true;
  }
  return to_return;
}

void delayMap::update_maps(const vector<bool>& slotmask) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < maps_by_slot.size(); i++) {
    if(slotmask[i]) {
      if(updated_to[i] < least_up_to_date) {
        least_up_to_date = updated_to[i];
      }
    }
  }
  if(current_site != least_up_to_date) {
    size_t suffixes_size = current_site - least_up_to_date;

    vector<DPUpdateMap> suffixes = 
              vector<DPUpdateMap>(suffixes_size, DPUpdateMap());
    suffixes[0] = maps_by_site[current_site];

    for(size_t i = 1; i < suffixes_size; i++) {
      suffixes[i] = suffixes[i-1].compose(maps_by_site[current_site - i]);
    }    
    
    for(size_t i = 0; i < maps_by_slot.size(); i++) {
      if(slotmask[i]) {
        if(updated_to[i] != current_site) {
          // j is the slot's index in the suffix-vector
          size_t j = current_site - updated_to[i] - 1;
          maps_by_slot[i] = suffixes[j].of(maps_by_slot[i]);
          updated_to[i] = current_site;
        }
      }
    }
  }
  updated_maps = true;
  return;
}

void delayMap::update_maps(const vector<size_t>& slots) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < slots.size(); i++) {
    if(updated_to[slots[i]] < least_up_to_date) {
      least_up_to_date = updated_to[slots[i]];
    }
  }
  if(current_site != least_up_to_date) {
    size_t suffixes_size = current_site - least_up_to_date;

    vector<DPUpdateMap> suffixes = 
              vector<DPUpdateMap>(suffixes_size, DPUpdateMap());
    suffixes[0] = maps_by_site[current_site];

    for(size_t i = 1; i < suffixes_size; i++) {
      suffixes[i] = suffixes[i-1].compose(maps_by_site[current_site - i]);
    }    
    
    for(size_t i = 0; i < slots.size(); i++) {
      if(updated_to[slots[i]] != current_site) {
        // j is the slot's index in the suffix-vector
        size_t j = current_site - updated_to[slots[i]] - 1;
        maps_by_slot[slots[i]] = suffixes[j].of(maps_by_slot[slots[i]]);
        updated_to[slots[i]] = current_site;
      }
    }
  }
  updated_maps = true;
  return;
}

void delayMap::increment_site_marker() {
  current_site++;
  added_span = false;
  updated_maps = false;
  return;
}

void delayMap::delete_slot(size_t slot) {
  maps_by_slot[slot] = DPUpdateMap(0);
  counts[slot] = 0;
  updated_to[slot] = current_site;
  empty_map_slots.push_back(slot);
  return;
}

void delayMap::decrement_slot(size_t slot) {
  if(counts[slot] == 1) {
    delete_slot(slot);
  } else {
    counts[slot]--;
  }
  return;
}

void delayMap::remove_row_from_slot(size_t row) {
  decrement_slot(slots_by_row[row]);
  // unassigned row is given max possible slot index + 1 to ensure that
  // accessing it will throw an error
  slots_by_row[row] = slots_by_row.size();
  return;
}

void delayMap::add_map(const DPUpdateMap& map) {
  if(empty_map_slots.size() == 0) {
    newest_index = maps_by_slot.size();
    maps_by_slot.push_back(map);
    counts.push_back(0);
    updated_to.push_back(current_site);
    return;
  } else {
    newest_index = empty_map_slots.back();
    empty_map_slots.pop_back();
    maps_by_slot[newest_index] = map;
    counts[newest_index] = 0;
    updated_to[newest_index] = current_site;
    return;
  }
}

void delayMap::add_map(double coefficient, double constant) {
  add_map(DPUpdateMap(coefficient, constant));
  return;
}

double delayMap::get_constant(size_t row) const {
  return maps_by_slot[slots_by_row[row]].constant;
}

double delayMap::get_coefficient(size_t row) const {
  return maps_by_slot[slots_by_row[row]].coefficient;
}

const DPUpdateMap& delayMap::get_map(size_t row) const {
  return maps_by_slot[slots_by_row[row]];
}

const vector<DPUpdateMap>& delayMap::get_maps() const {
  return maps_by_slot;
}

vector<DPUpdateMap>& delayMap::get_maps() {
  return maps_by_slot;
}

const vector<size_t>& delayMap::get_map_indices() const {
  return slots_by_row;
}

void delayMap::update_map_with_span(const DPUpdateMap& span_map) {
  add_map_for_site(span_map);
  return;
}

void delayMap::update_map_with_span(double coefficient, double constant) {
  update_map_with_span(DPUpdateMap(coefficient, constant));
  return;
}

void delayMap::add_map_for_site(const DPUpdateMap& site_map) {
  increment_site_marker();
  maps_by_site.push_back(site_map);
  return;
}

void delayMap::add_map_for_site(double coefficient, double constant) {
  add_map_for_site(DPUpdateMap(coefficient, constant));
  return;
}

size_t delayMap::last_update(size_t row) {
  if(slots_by_row[row] != slots_by_row.size()) {
    return updated_to[slots_by_row[row]];
  } else {
    return current_site;
  }
}

const vector<DPUpdateMap>& delayMap::get_maps_by_site() const {
  return maps_by_site.get_elements();
}

void delayMap::reset_rows(const vector<size_t>& rows) {
  for(size_t i = 0; i < rows.size(); i++) {
    remove_row_from_slot(rows[i]);
  }
  add_identity_map();
  for(size_t i = 0; i < rows.size(); i++) {
    assign_row_to_newest_index(rows[i]);
  }
}

void delayMap::reset_rows(const rowSet& rows) {
  for(size_t i = 0; i < rows.size(); i++) {
    remove_row_from_slot(rows[i]);
  }
  add_identity_map();
  for(size_t i = 0; i < rows.size(); i++) {
    assign_row_to_newest_index(rows[i]);
  }
}

void delayMap::update_map_with_active_rows(const vector<size_t>& active_rows) {
  vector<size_t> slots = rows_to_slots(active_rows);
  update_maps(slots);
}

void delayMap::update_map_with_active_rows(const rowSet& active_rows) {
  update_maps(rows_to_slotmask(active_rows));
}

size_t delayMap::number_of_slots() const {
  return counts.size() - empty_map_slots.size();
}

size_t delayMap::row_updated_to(size_t row) const {
  return updated_to[slots_by_row[row]];
}

size_t delayMap::get_current_site() const {
  return current_site;
}

size_t delayMap::get_slot(size_t row) const {
  return slots_by_row[row];
}

double delayMap::evaluate(size_t row, double value) const {
  return maps_by_slot[slots_by_row[row]].of(value);
}
