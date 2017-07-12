#include "lh_delay_multiplier.hpp"
#include "lh_math.hpp"

using namespace std;

delayMap::delayMap(size_t rows, size_t start) : 
            dM_start(start) {
  slots_by_row = vector<size_t>(rows, rows);
  dM_start = start;
  current_site = start;
}

void delayMap::assign_row_to_newest_index(size_t row) {
  //TODO: complain if slots_by_row[row] != |H|
  slots_by_row[row] = newest_index;
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

void delayMap::update_maps(const vector<size_t>& indices) {
  size_t oldest_start_position = current_site;
  for(size_t i = 0; i < indices.size(); i++) {
    if(updated_to[indices[i]] < oldest_start_position) {
      oldest_start_position = updated_to[indices[i]];
    }
  }
  if(current_site > oldest_start_position) {
    vector<DPUpdateMap> suffixes;
    suffixes.push_back(maps_by_site[current_site]);
    for(size_t i = 1; i <= current_site - oldest_start_position; i++) {
      suffixes.push_back(
                  suffixes.back().compose(maps_by_site[current_site - i]));
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
  maps_by_slot[slot] = DPUpdateMap(0,0);
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

void delayMap::add_map(DPUpdateMap map) {
  if(empty_map_slots.size() == 0) {
    newest_index = maps_by_slot.size();
    maps_by_slot.push_back(map);
    counts.push_back(1);
    updated_to.push_back(current_site);
    return;
  } else {
    newest_index = empty_map_slots.back();
    empty_map_slots.pop_back();
    maps_by_slot[newest_index] = map;
    counts[newest_index] = 1;
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

vector<DPUpdateMap> delayMap::get_maps() {
  return maps_by_slot;
}

vector<size_t> delayMap::get_map_indices() {
  return slots_by_row;
}

void delayMap::update_map_with_span(DPUpdateMap span_map) {
  added_span = true;
  maps_by_site.back() = span_map.compose(maps_by_site.back());
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