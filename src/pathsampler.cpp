#include "pathsampler.hpp"
#include <random>
#include <unordered_set>
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace std;

PathSampler::PathSampler(size_t total_number) 
	: total_number(total_number)
{}

void PathSampler::select_single_subset(vector<size_t>& result, size_t sample_size) const {
	assert(sample_size <= this->total_number);
	unordered_set<int> sample;
	default_random_engine generator;

	for(size_t d = this->total_number - sample_size; d < this->total_number; ++d) {
		int t = uniform_int_distribution<>(0, d)(generator);
		if (sample.find(t) == sample.end() )
			sample.insert(t);
		else
			sample.insert(d);
	}
	std::copy(sample.begin(), sample.end(), std::back_inserter(result));
	sort(result.begin(), result.end());
}

void PathSampler::select_multiple_subsets(vector<vector<size_t>>& result, size_t sample_size, size_t n) const {
	for (size_t i = 0; i < n; ++i) {
		vector<size_t> sample;
		this->select_single_subset(sample, sample_size);
		result.push_back(sample);
	}
}

void PathSampler::partition_paths(vector<vector<size_t>>& result, size_t sample_size) const {
	vector<size_t> all_paths = {};
	for (size_t i = 0; i < this->total_number; ++i) {
		all_paths.push_back(i);
	}
	random_shuffle(all_paths.begin(), all_paths.end());
	
	// partition into parts of length sample_size
	for (size_t i = 0; i < all_paths.size(); i += sample_size) {
		auto last = min(all_paths.size(), i+sample_size);
		vector<size_t> subset(all_paths.begin()+i, all_paths.begin()+last);
		sort(subset.begin(), subset.end());
		result.push_back(subset);
	}

	// if the last sample is shorter than the sample_size, randomly add additional paths
	size_t missing = sample_size - result.at(result.size()-1).size();
	if (missing > 0) {
		this->select_single_subset(result[result.size()-1], missing);
	}
}
