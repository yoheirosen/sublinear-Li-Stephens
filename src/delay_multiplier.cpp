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

lazyEvalMap::lazyEvalMap() {
  
}

void lazyEvalMap::increment_site_marker() {
  current_site++;
}

lazyEvalMap::lazyEvalMap(size_t rows, size_t start) : 
	row_to_mapclass(vector<size_t>(rows, 0)), 
	mapclass_size(vector<size_t>(1, rows)),
	current_site(start),
	mapclass_to_last_updated(vector<size_t>(1, start)),
	newest_mapclass(0),
	mapclass_to_map(vector<DPUpdateMap>(1, DPUpdateMap(0))),
	maps_by_site(mapHistory(DPUpdateMap(0), start)) {

}

void lazyEvalMap::add_identity_map() {
  add_map(DPUpdateMap(0));
  return;
}

lazyEvalMap::lazyEvalMap(const lazyEvalMap &other) {
	current_site = other.current_site;
	row_to_mapclass = other.row_to_mapclass;
	mapclass_to_last_updated = other.mapclass_to_last_updated;
  size_t oldest_seen = current_site;
  for(size_t i = 0; i < mapclass_to_last_updated.size(); i++) {
    if(mapclass_to_last_updated[i] < oldest_seen) {
      oldest_seen = mapclass_to_last_updated[i];
    }
  }
  maps_by_site = mapHistory(other.maps_by_site, oldest_seen);
	mapclass_to_map = other.mapclass_to_map;
	mapclass_size = other.mapclass_size;
	empty_mapclass_indices = other.empty_mapclass_indices;
}

void lazyEvalMap::assign_row_to_newest_mapclass(size_t row) {
  //TODO: complain if row_to_mapclass[row] != |H|
  row_to_mapclass[row] = newest_mapclass;
  mapclass_size[newest_mapclass]++;
  return;
}

void lazyEvalMap::hard_clear_all() {
  for(int i = 0; i < mapclass_to_map.size(); i++) {
    delete_mapclass(i);
  }
  add_identity_map();
  for(int i = 0; i < row_to_mapclass.size(); i++) {
    assign_row_to_newest_mapclass(i);
  }
  return;
}

void lazyEvalMap::hard_update_all() {
  vector<size_t> non_empty_mapclasses;
  
  for(int i = 0; i < mapclass_size.size(); i++) {
    if(mapclass_size[i] != 0) {
      non_empty_mapclasses.push_back(i);
    }
  }
  update_maps(non_empty_mapclasses);
  return;
}

vector<size_t> lazyEvalMap::rows_to_mapclasses(const rowSet& rows) const {
  size_t rows_size = rows.size();
  vector<char> seen = vector<char>(mapclass_to_map.size(), 0);
  vector<size_t> to_return;
  rowSet::const_iterator it = rows.begin();
  for(it; it != rows.end(); ++it) {
    if(seen[row_to_mapclass[*it]] == 0) {
      to_return.push_back(row_to_mapclass[*it]);
    }
    seen[row_to_mapclass[*it]] = 1;
  }
  return to_return;
}

vector<bool> lazyEvalMap::rows_to_mapclassmask(const rowSet& rows) const {
  size_t rows_size = rows.size();
  vector<bool> to_return = vector<bool>(mapclass_to_map.size(), false);
  for(int i = 0; i < rows_size; i++) {
    to_return[row_to_mapclass[rows[i]]] = true;
  }
  return to_return;
}

void lazyEvalMap::update_maps(const vector<bool>& mapclassmask) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < mapclass_to_map.size(); i++) {
    if(mapclassmask[i]) {
      if(mapclass_to_last_updated[i] < least_up_to_date) {
        least_up_to_date = mapclass_to_last_updated[i];
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
    
    for(size_t i = 0; i < mapclass_to_map.size(); i++) {
      if(mapclassmask[i]) {
        if(mapclass_to_last_updated[i] != current_site) {
          // j is the mapclass's index in the suffix-vector
          size_t j = current_site - mapclass_to_last_updated[i] - 1;
          mapclass_to_map[i] = suffixes[j].of(mapclass_to_map[i]);
          mapclass_to_last_updated[i] = current_site;
        }
      }
    }
  }
  return;
}

void lazyEvalMap::update_maps(const vector<size_t>& mapclasses) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < mapclasses.size(); i++) {
    if(mapclass_to_last_updated[mapclasses[i]] < least_up_to_date) {
      least_up_to_date = mapclass_to_last_updated[mapclasses[i]];
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
    
    for(size_t i = 0; i < mapclasses.size(); i++) {
      if(mapclass_to_last_updated[mapclasses[i]] != current_site) {
        // j is the mapclass's index in the suffix-vector
        size_t j = current_site - mapclass_to_last_updated[mapclasses[i]] - 1;
        mapclass_to_map[mapclasses[i]] = suffixes[j].of(mapclass_to_map[mapclasses[i]]);
        mapclass_to_last_updated[mapclasses[i]] = current_site;
      }
    }
  }
  return;
}

void lazyEvalMap::delete_mapclass(size_t mapclass) {
  // mapclass_to_map[mapclass] = DPUpdateMap(0);
  mapclass_size[mapclass] = 0;
  mapclass_to_last_updated[mapclass] = current_site;
  empty_mapclass_indices.push_back(mapclass);
  return;
}

void lazyEvalMap::decrement_mapclass(size_t mapclass) {
  if(mapclass_size[mapclass] == 1) {
    delete_mapclass(mapclass);
  } else {
    --mapclass_size[mapclass];
  }
  return;
}

void lazyEvalMap::remove_row_from_mapclass(size_t row) {
  decrement_mapclass(row_to_mapclass[row]);
  // unassigned row is given max possible mapclass index + 1 to ensure that
  // accessing it will throw an error
  row_to_mapclass[row] = row_to_mapclass.size();
  return;
}

void lazyEvalMap::add_map(const DPUpdateMap& map) {
  if(empty_mapclass_indices.size() == 0) {
    newest_mapclass = mapclass_to_map.size();
    mapclass_to_map.push_back(map);
    mapclass_size.push_back(0);
    mapclass_to_last_updated.push_back(current_site);
    return;
  } else {
    newest_mapclass = empty_mapclass_indices.back();
    empty_mapclass_indices.pop_back();
    mapclass_to_map[newest_mapclass] = map;
    mapclass_size[newest_mapclass] = 0;
    mapclass_to_last_updated[newest_mapclass] = current_site;
    return;
  }
}

double lazyEvalMap::get_constant(size_t row) const {
  return mapclass_to_map[row_to_mapclass[row]].constant;
}

double lazyEvalMap::get_coefficient(size_t row) const {
  return mapclass_to_map[row_to_mapclass[row]].coefficient;
}

const DPUpdateMap& lazyEvalMap::get_map(size_t row) const {
  return mapclass_to_map[row_to_mapclass[row]];
}

const vector<DPUpdateMap>& lazyEvalMap::get_maps() const {
  return mapclass_to_map;
}

vector<DPUpdateMap>& lazyEvalMap::get_maps() {
  return mapclass_to_map;
}

const vector<size_t>& lazyEvalMap::get_map_indices() const {
  return row_to_mapclass;
}

void lazyEvalMap::update_map_with_span(const DPUpdateMap& span_map) {
  add_map_for_site(span_map);
  return;
}

void lazyEvalMap::add_map_for_site(const DPUpdateMap& site_map) {
  current_site++;
  maps_by_site.push_back(site_map);
  return;
}

size_t lazyEvalMap::last_update(size_t row) const {
  if(row_to_mapclass[row] != row_to_mapclass.size()) {
    return mapclass_to_last_updated[row_to_mapclass[row]];
  } else {
    return current_site;
  }
}

const vector<DPUpdateMap>& lazyEvalMap::get_maps_by_site() const {
  return maps_by_site.get_elements();
}

void lazyEvalMap::reset_rows(const rowSet& rows) {
  size_t rows_size = rows.size();
  for(size_t i = 0; i < rows_size; i++) {
    remove_row_from_mapclass(rows[i]);
  }
  add_identity_map();
  for(size_t i = 0; i < rows_size; i++) {
    assign_row_to_newest_mapclass(rows[i]);
  }
}

void lazyEvalMap::update_map_with_active_rows(const rowSet& active_rows) {
  update_maps(rows_to_mapclassmask(active_rows));
}

size_t lazyEvalMap::number_of_mapclasses() const {
  return mapclass_size.size() - empty_mapclass_indices.size();
}

size_t lazyEvalMap::row_updated_to(size_t row) const {
  return mapclass_to_last_updated[row_to_mapclass[row]];
}

size_t lazyEvalMap::get_current_site() const {
  return current_site;
}

size_t lazyEvalMap::get_mapclass(size_t row) const {
  return row_to_mapclass[row];
}

double lazyEvalMap::evaluate(size_t row, double value) const {
  return mapclass_to_map[row_to_mapclass[row]].of(value);
}
