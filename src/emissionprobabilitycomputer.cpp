#include <stdexcept>
#include <sstream>
#include "emissionprobabilitycomputer.hpp"

using namespace std;

pair<size_t,size_t> sorted_paths(size_t p1, size_t p2) {
	if (p1 < p2) {
		return make_pair(p1, p2);
	} else {
		return make_pair(p2, p1);
	}
}

EmissionProbabilityComputer::EmissionProbabilityComputer(UniqueKmers* uniquekmers)
	:uniquekmers(uniquekmers)
{
	// precompute all probabilities
	vector<size_t> path_ids;
	this->uniquekmers->get_path_ids(path_ids);
	size_t nr_paths = path_ids.size();
	for (size_t i = 0; i < nr_paths; ++i) {
		for (size_t j = i; j < nr_paths; ++j) {
			size_t path_id1 = path_ids[i];
			size_t path_id2 = path_ids[j];
			long double emission_prob = compute_emission_probability(path_id1, path_id2);

			// store probability (sort path ids numerically, to store each combination only once
			pair<size_t,size_t> paths = sorted_paths(path_id1, path_id2);
			this->paths_to_prob[paths] = emission_prob;
		}
	}
}

long double EmissionProbabilityComputer::get_emission_probability(int path_id1, int path_id2) const {
	pair<size_t,size_t> paths = sorted_paths(path_id1, path_id2);
	auto it = this->paths_to_prob.find(paths);
	if (it != this->paths_to_prob.end()) {
		return this->paths_to_prob.at(paths);
	} else {
		throw runtime_error("EmissionProbabilityComputer::get_emission_probability: invalid path ids given.");
	}
}

long double EmissionProbabilityComputer::compute_emission_probability(int path_id1, int path_id2){
	long double result = 1.0L;
	// combine the two paths to get expected kmer copy numbers
	CopyNumberAssignment cna = this->uniquekmers->combine_paths(path_id1, path_id2);
	for (size_t i = 0; i < this->uniquekmers->size(); ++i){
		unsigned int expected_kmer_count = cna.get_position(i);
		// multiply result by probability of expected kmer count
		result *= uniquekmers->kmer_to_copynumber[i].get_probability_of(expected_kmer_count);
	}
	return result;
}
