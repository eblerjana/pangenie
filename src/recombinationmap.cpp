#include <stdexcept>
#include <utility>
#include <sstream>
#include "recombinationmap.hpp"

using namespace std;

float interpolate (size_t point, size_t start_pos, size_t end_pos, float start_val, float end_val) {
	assert (start_pos <= point);
	assert (point <= end_pos);

	if ( (start_pos == point) && (point == end_pos) ) return start_val;
	return start_val + ((point - start_pos) * (end_val - start_val) / (end_pos - start_pos));
}

void load_genetic_map(string filename, vector<pair<size_t,float>>* result) {
	ifstream file(filename);
	if (!file.good()) {
		throw runtime_error("RecombinationMap::load_genetic_map: file " + filename + " cannot be opened.");
	}

	string line;
	while(getline(file, line)) {
		// check if header line
		if (line.substr(0,8) == "position") continue;

		// parse the fields
		vector<string> fields;

		if (line == "") continue;

		string token;
		istringstream iss (line);
		while (getline(iss, token, ' ')) {
			fields.push_back(token);
		}

		// make sure the file is properly formatted
		if (fields.size() != 3) {
			throw runtime_error("RecombinationMap::load_genetic_map: genetic map file is malformed.");
		}

		// convert strings to numbers
		stringstream ss_pos(fields[0]);
		size_t position;
		ss_pos >> position;

		stringstream ss_dist(fields[2]);
		float dist;
		ss_dist >> dist;

		// store position and cummulative distance
		result->push_back(make_pair(position, dist));
	}
}

// TODO: this needs access to the variant positions also!!
RecombinationMap::RecombinationMap(std::string& filename, VariantReader* variants, string chromosome)
	:uniform(false)
{
	// read the recombination map file
	vector<pair<size_t, float>> recomb_map;
	load_genetic_map(filename, &recomb_map);
	compute_recombination_cost_map(&recomb_map, variants, chromosome);

}

RecombinationMap::RecombinationMap(float recomb_rate)
	: uniform(true),
	 recomb_rate((long double) recomb_rate)

{}


void RecombinationMap::compute_recombination_cost_map(vector<pair<size_t,float>>* genetic_map, VariantReader* variants, string chromosome) {
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
		float d = 0.0;
		if (i_none) {
			assert (!j_none);
			d = interpolate(position, 0, genetic_map->at(j).first, 0,  genetic_map->at(j).second);
		} else if (j_none) {
			// point is outside of the genetic map
			float avg_rate = genetic_map->at(map_size-1).second / genetic_map->at(map_size -1).first;
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
long double RecombinationMap::compute_recombination_probability(size_t left_variant, size_t left_variant_id, size_t right_variant, size_t right_variant_id) {

	if (right_variant != (left_variant + 1)) {
		throw runtime_error("RecombinationMap::compute_recombination_probability: variant ids must be consecutive.");
	}

	if (right_variant >= this->cumulative_distances.size()) {
		throw runtime_error("RecombinationMap::compute_recombination_probability: right variant id exeeds the number of variants stored.");
	}

	if (this->uniform) {
		return (right_variant - left_variant) * 0.000001 * (this->recomb_rate);
	} else {
		if (left_variant == 0) return 0.0L;
		return this->cumulative_distances[right_variant_id] - cumulative_distances[left_variant_id];
	}
}



