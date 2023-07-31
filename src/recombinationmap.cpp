#include <stdexcept>
#include <utility>
#include <sstream>
#include "recombinationmap.hpp"

using namespace std;

long double interpolate (size_t point, size_t start_pos, size_t end_pos, long double start_val, long double end_val) {
	assert (start_pos <= point);
	assert (point <= end_pos);

	if ( (start_pos == point) && (point == end_pos) ) return start_val;
	return start_val + ((point - start_pos) * (end_val - start_val) / (end_pos - start_pos));
}

void parse_map_inputs(string filename, map<string, string>& output_maps) {
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("parse_map_inputs: file " + filename + " cannot be opened.");
	}

	string line;
	while(getline(file, line)){
		vector<string> fields;
		if (line == "") continue;
		string token;
		istringstream iss (line);
		while (getline(iss, token, '\t')){
			fields.push_back(token);
		}

		if (fields.size() != 2) {
			throw runtime_error("parse_map_inputs: map config file is malformed.");
		}

		output_maps[fields[0]] = fields[1];
	}
}

void load_genetic_map(string filename, vector<pair<size_t,long double>>* result) {
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("load_genetic_map: file " + filename + " cannot be opened.");
	}

	string line;
	while(getline(file, line)) {
		// check if header line
		if (line.substr(0,3) == "pos") continue;

		// parse the fields
		vector<string> fields;
		if (line == "") continue;
		string token;
		istringstream iss (line);
		while (getline(iss, token, '\t')) {
			fields.push_back(token);
		}

		// make sure the file is properly formatted
		if (fields.size() != 3) {
			throw runtime_error("load_genetic_map: genetic map file is malformed.");
		}

		// convert strings to numbers
		stringstream ss_pos(fields[0]);
		size_t position;
		ss_pos >> position;

		stringstream ss_dist(fields[2]);
		long double dist;
		ss_dist >> dist;

		// store position and cummulative distance
		result->push_back(make_pair(position, dist));
	}
}

RecombinationMap::RecombinationMap(std::string& filename, VariantReader* variants, string chromosome)
	:uniform(false)
{
	// read the recombination map file
	vector<pair<size_t, long double>> recomb_map;
	load_genetic_map(filename, &recomb_map);
	compute_recombination_cost_map(&recomb_map, variants, chromosome);

}

RecombinationMap::RecombinationMap(long double recomb_rate)
	: uniform(true),
	 recomb_rate((long double) recomb_rate)

{}


void RecombinationMap::compute_recombination_cost_map(vector<pair<size_t,long double>>* genetic_map, VariantReader* variants, string chromosome) {
	// NOTE: c++ implementation of the WhatsHap code: https://github.com/whatshap/whatshap/blob/main/whatshap/pedigree.py
	assert (genetic_map->size() > 0);

	// i and j are such that genetic_map[i].position <= position <= genetic_map[j].position
	// i and j are None if no such values exist (because we are at the end of the list)
	size_t i = 0;
	size_t j = 0;
	bool i_none = true;
	bool j_none = false;

	size_t nr_variants = variants->size_of(chromosome);
	size_t map_size = genetic_map->size();

	for (size_t v = 0; v < nr_variants; ++v) {
		size_t position = variants->get_variant(chromosome, v).get_start_position();

		if (i_none && (genetic_map->at(0).first <= position)) {
			i = 0;
			i_none = false;
		}

		while (!i_none && ( (i+1) < map_size ) && (genetic_map->at(i+1).first <= position) ) {
			i += 1;
		}

		// update j to meet the invariant
		while ( (!j_none) && (genetic_map->at(j).first < position) ) {
			if (j+1 < map_size) {
				j += 1;
			} else {
				j_none = true;
			}
		}

		// interpolate
		long double d = 0.0;
		if (i_none) {
			assert (!j_none);
			d = interpolate(position, 0, genetic_map->at(j).first, 0,  genetic_map->at(j).second);
		} else if (j_none) {
			// point is outside of the genetic map
			long double avg_rate = genetic_map->at(map_size-1).second / genetic_map->at(map_size -1).first;
			d = genetic_map->at(map_size-1).second + (position - genetic_map->at(map_size-1).first) * avg_rate;
		} else {
			assert (genetic_map->at(i).first <= position);
			assert (position <= genetic_map->at(j).first);
			d = interpolate (position, genetic_map->at(i).first, genetic_map->at(j).first, genetic_map->at(i).second, genetic_map->at(j).second);
		}

		this->cumulative_distances.push_back(d);
	}
}


// TODO: are IDs available?? --> Yes, should be.
long double RecombinationMap::compute_recombination_probability(size_t left_variant, size_t left_variant_id, size_t right_variant, size_t right_variant_id) const {

	if (left_variant_id >= right_variant_id) {
		throw runtime_error("RecombinationMap::compute_recombination_probability: left_variant_id must be smaller than right_variant_id.");
	}

	if (left_variant > right_variant) {
		throw runtime_error("RecombinationMap::compute_recombination_probability: left variant position is larger than right variant position.");
	}

	if (this->uniform) {
		return (right_variant - left_variant) * 0.000001 * (this->recomb_rate);
	} else {
		long double minimum_genetic_distance = 1e-10L;
		if (right_variant_id >= this->cumulative_distances.size()) {
			throw runtime_error("RecombinationMap::compute_recombination_probability: right variant id exeeds the number of variants stored.");
		}
	//	if (left_variant_id == 0) return 0.0L;
		return max(this->cumulative_distances[right_variant_id] - cumulative_distances[left_variant_id], minimum_genetic_distance);
	}
}



