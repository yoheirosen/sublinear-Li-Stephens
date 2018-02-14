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
	row_to_eqclass(vector<size_t>(rows, 0)), 
	eqclass_size(vector<size_t>(1, rows)),
	current_site(start),
	eqclass_last_updated(vector<size_t>(1, start)),
	newest_eqclass(0),
	eqclass_to_map(vector<DPUpdateMap>(1, DPUpdateMap(0))),
	map_history(mapHistory(DPUpdateMap(0), start)) {

}

void lazyEvalMap::add_identity_eqclass() {
  add_eqclass(DPUpdateMap(0));
  return;
}

lazyEvalMap::lazyEvalMap(const lazyEvalMap &other) {
	current_site = other.current_site;
	row_to_eqclass = other.row_to_eqclass;
	eqclass_last_updated = other.eqclass_last_updated;
  size_t oldest_seen = current_site;
  for(size_t i = 0; i < eqclass_last_updated.size(); i++) {
    if(eqclass_last_updated[i] < oldest_seen) {
      oldest_seen = eqclass_last_updated[i];
    }
  }
  map_history = mapHistory(other.map_history, oldest_seen);
	eqclass_to_map = other.eqclass_to_map;
	eqclass_size = other.eqclass_size;
	empty_eqclass_indices = other.empty_eqclass_indices;
}

void lazyEvalMap::assign_row_to_newest_eqclass(size_t row) {
  //TODO: complain if row_to_eqclass[row] != |H|
  row_to_eqclass[row] = newest_eqclass;
  eqclass_size[newest_eqclass]++;
  return;
}

void lazyEvalMap::hard_clear_all() {
  for(int i = 0; i < eqclass_to_map.size(); i++) {
    delete_eqclass(i);
  }
  add_identity_eqclass();
  for(int i = 0; i < row_to_eqclass.size(); i++) {
    assign_row_to_newest_eqclass(i);
  }
  return;
}

void lazyEvalMap::hard_update_all() {
  vector<size_t> non_empty_eqclasses;
  
  for(int i = 0; i < eqclass_size.size(); i++) {
    if(eqclass_size[i] != 0) {
      non_empty_eqclasses.push_back(i);
    }
  }
  update_maps(non_empty_eqclasses);
  return;
}

vector<size_t> lazyEvalMap::rows_to_eqclasses(const rowSet& rows) const {
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

void lazyEvalMap::update_maps(const vector<size_t>& eqclasses) {
  size_t least_up_to_date = current_site;
  for(size_t i = 0; i < eqclasses.size(); i++) {
    if(eqclass_last_updated[eqclasses[i]] < least_up_to_date) {
      least_up_to_date = eqclass_last_updated[eqclasses[i]];
    }
  }
  if(current_site != least_up_to_date) {
    size_t suffixes_size = current_site - least_up_to_date;

    vector<DPUpdateMap> suffixes = 
              vector<DPUpdateMap>(suffixes_size, DPUpdateMap());
    suffixes[0] = map_history[current_site];

    for(size_t i = 1; i < suffixes_size; i++) {
      suffixes[i] = suffixes[i-1].compose(map_history[current_site - i]);
    }    
    
    for(size_t i = 0; i < eqclasses.size(); i++) {
      if(eqclass_last_updated[eqclasses[i]] != current_site) {
        // j is the eqclass's index in the suffix-vector
        size_t j = current_site - eqclass_last_updated[eqclasses[i]] - 1;
        eqclass_to_map[eqclasses[i]] = suffixes[j].of(eqclass_to_map[eqclasses[i]]);
        eqclass_last_updated[eqclasses[i]] = current_site;
      }
    }
  }
  return;
}

void lazyEvalMap::delete_eqclass(size_t eqclass) {
  // eqclass_to_map[eqclass] = DPUpdateMap(0);
  eqclass_size[eqclass] = 0;
  eqclass_last_updated[eqclass] = current_site;
  empty_eqclass_indices.push_back(eqclass);
  return;
}

void lazyEvalMap::decrement_eqclass(size_t eqclass) {
  if(eqclass_size[eqclass] == 1) {
    delete_eqclass(eqclass);
  } else {
    --eqclass_size[eqclass];
  }
  return;
}

void lazyEvalMap::remove_row_from_eqclass(size_t row) {
  decrement_eqclass(row_to_eqclass[row]);
  // unassigned row is given max possible eqclass index + 1 to ensure that
  // accessing it will throw an error
  row_to_eqclass[row] = row_to_eqclass.size();
  return;
}

void lazyEvalMap::add_eqclass(const DPUpdateMap& map) {
  if(empty_eqclass_indices.size() == 0) {
    newest_eqclass = eqclass_to_map.size();
    eqclass_to_map.push_back(map);
    eqclass_size.push_back(0);
    eqclass_last_updated.push_back(current_site);
    return;
  } else {
    newest_eqclass = empty_eqclass_indices.back();
    empty_eqclass_indices.pop_back();
    eqclass_to_map[newest_eqclass] = map;
    eqclass_size[newest_eqclass] = 0;
    eqclass_last_updated[newest_eqclass] = current_site;
    return;
  }
}

double lazyEvalMap::get_constant(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]].constant;
}

double lazyEvalMap::get_coefficient(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]].coefficient;
}

const DPUpdateMap& lazyEvalMap::get_map(size_t row) const {
  return eqclass_to_map[row_to_eqclass[row]];
}

const vector<DPUpdateMap>& lazyEvalMap::get_maps() const {
  return eqclass_to_map;
}

vector<DPUpdateMap>& lazyEvalMap::get_maps() {
  return eqclass_to_map;
}

const vector<size_t>& lazyEvalMap::get_map_indices() const {
  return row_to_eqclass;
}

void lazyEvalMap::stage_map_for_span(const DPUpdateMap& span_map) {
  stage_map_for_site(span_map);
  return;
}

void lazyEvalMap::stage_map_for_site(const DPUpdateMap& site_map) {
  current_site++;
  map_history.push_back(site_map);
  return;
}

size_t lazyEvalMap::last_update(size_t row) const {
  if(row_to_eqclass[row] != row_to_eqclass.size()) {
    return eqclass_last_updated[row_to_eqclass[row]];
  } else {
    return current_site;
  }
}

const vector<DPUpdateMap>& lazyEvalMap::get_map_history() const {
  return map_history.get_elements();
}
void lazyEvalMap::reset_rows(const rowSet& rows) {
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

void lazyEvalMap::update_active_rows(const rowSet& active_rows) {
  // update_maps(rows_to_eqclassmask(active_rows));
  update_maps(rows_to_eqclasses(active_rows));
}

size_t lazyEvalMap::number_of_eqclasses() const {
  return eqclass_size.size() - empty_eqclass_indices.size();
}

size_t lazyEvalMap::row_updated_to(size_t row) const {
  return eqclass_last_updated[row_to_eqclass[row]];
}

size_t lazyEvalMap::get_current_site() const {
  return current_site;
}

size_t lazyEvalMap::get_eqclass(size_t row) const {
  return row_to_eqclass[row];
}

double lazyEvalMap::evaluate(size_t row, double value) const {
  return eqclass_to_map[row_to_eqclass[row]].of(value);
}
