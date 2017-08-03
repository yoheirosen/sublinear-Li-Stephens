#include "lh_delay_multiplier.hpp"
#include "lh_math.hpp"

using namespace std;

delayMap::delayMap() {
  
}

delayMap::delayMap(size_t rows, size_t start) : 
            dM_start(start) {
  current_site = start;
  add_identity_map();
  slots_by_row = vector<size_t>(rows, 0);
  counts = {rows};
  dM_start = start;
}

delayMap::delayMap(const delayMap &other) {
	added_span = other.added_span;
	updated_maps = other.updated_maps;
	dM_start = other.dM_start;
	current_site = other.current_site;
	maps_by_site = other.maps_by_site;
	slots_by_row = other.slots_by_row;
	updated_to = other.updated_to;
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
  vector<size_t> indices;
  for(int i = 0; i < counts.size(); i++) {
    if(counts[i] != 0) {
      indices.push_back(i);
    }
  }
  update_maps(indices);
  return;
}

vector<size_t> delayMap::rows_to_slots(const vector<size_t>& rows) {
  vector<bool> seen = vector<bool>(maps_by_site.size(), false);
  vector<size_t> to_return;
  for(int i = 0; i < rows.size(); i++) {
    if(!(seen[slots_by_row[rows[i]]])) {
      to_return.push_back(slots_by_row[rows[i]]);
    }
    seen[slots_by_row[rows[i]]] = true;
  }
  return to_return;
}

void delayMap::update_maps(const vector<size_t>& indices) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < indices.size(); i++) {
    if(updated_to[indices[i]] < least_up_to_date) {
      least_up_to_date = updated_to[indices[i]];
    }
  }
  if(current_site > least_up_to_date) {
    vector<DPUpdateMap> suffixes;
    
    // there is no map for the first site
    suffixes.push_back(maps_by_site[current_site - 1]);

    size_t number_of_suffixes = current_site - least_up_to_date;
    for(size_t i = 1; i < number_of_suffixes; i++) {
      suffixes.push_back(
                  suffixes.back().compose(maps_by_site[current_site - i - 1]));
    }
    
    for(size_t i = 0; i < indices.size(); i++) {
      size_t suffix_index = current_site - updated_to[indices[i]] - 1;
      maps_by_slot[indices[i]] = 
                suffixes[suffix_index].compose(maps_by_slot[indices[i]]);
      updated_to[indices[i]] = current_site;
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
  maps_by_slot[slot] = DPUpdateMap(0,-1);
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

void delayMap::add_identity_map() {
  add_map(0, -1);
  return;
}

void delayMap::add_map(DPUpdateMap map) {
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

double delayMap::get_constant(size_t row) {
  return maps_by_slot[slots_by_row[row]].constant;
}

double delayMap::get_coefficient(size_t row) {
  return maps_by_slot[slots_by_row[row]].coefficient;
}

DPUpdateMap delayMap::get_map(size_t row) {
  return maps_by_slot[slots_by_row[row]];
}


vector<DPUpdateMap> delayMap::get_maps() {
  return maps_by_slot;
}

vector<size_t> delayMap::get_map_indices() {
  return slots_by_row;
}

void delayMap::update_map_with_span(DPUpdateMap span_map) {
  add_map_for_site(span_map);
  return;
}

void delayMap::update_map_with_span(double coefficient, double constant) {
  update_map_with_span(DPUpdateMap(coefficient, constant));
  return;
}

void delayMap::add_map_for_site(DPUpdateMap site_map) {
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

vector<DPUpdateMap> delayMap::get_maps_by_site() {
  return maps_by_site;
}

void delayMap::reset_rows(vector<size_t>& rows) {
  for(size_t i = 0; i < rows.size(); i++) {
    remove_row_from_slot(rows[i]);
  }
  add_identity_map();
  for(size_t i = 0; i < rows.size(); i++) {
    assign_row_to_newest_index(rows[i]);
  }
}

void delayMap::update_map_with_active_rows(vector<size_t>& active_rows) {
  vector<size_t> slots = rows_to_slots(active_rows);
  update_maps(slots);
}