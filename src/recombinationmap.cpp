#include <stdexcept>
#include <utility>
#include <sstream>
#include "recombinationmap.hpp"

using namespace std;

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
RecombinationMap::RecombinationMap(std::string& filename)
	:uniform(false)
{
	// TODO: read the recombination map file and compute recombination probabilities
}

RecombinationMap::RecombinationMap(float recomb_rate)
	:recomb_rate((long double) recomb_rate),
	 uniform(true)
{}

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
		return this->cumulative_distances[right_variant_id] - cumulative_distances[left_variant_id];
	}
}



