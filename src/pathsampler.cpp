#include "pathsampler.hpp"
#include <random>
#include <unordered_set>
#include <iostream>
#include <cassert>
#include <algorithm>

using namespace std;

PathSampler::PathSampler(unsigned short total_number) 
	: total_number(total_number)
{}

void PathSampler::select_single_subset(vector<unsigned short>& result, unsigned short sample_size) const {
	assert(sample_size <= this->total_number);
	unordered_set<int> sample;
	default_random_engine generator;

	for(unsigned short d = this->total_number - sample_size; d < this->total_number; ++d) {
		int t = uniform_int_distribution<>(0, d)(generator);
		if (sample.find(t) == sample.end() )
			sample.insert(t);
		else
			sample.insert(d);
	}
	std::copy(sample.begin(), sample.end(), std::back_inserter(result));
	sort(result.begin(), result.end());
}

void PathSampler::select_multiple_subsets(vector<vector<unsigned short>>& result, unsigned short sample_size,  unsigned short n) const {
	for (unsigned short i = 0; i < n; ++i) {
		vector<unsigned short> sample;
		this->select_single_subset(sample, sample_size);
		result.push_back(sample);
	}
}

void PathSampler::partition_paths(vector<vector<unsigned short>>& result, unsigned short sample_size) const {
	vector<unsigned short> all_paths = {};
	for (unsigned short i = 0; i < this->total_number; ++i) {
		all_paths.push_back(i);
	}
	random_shuffle(all_paths.begin(), all_paths.end());
	
	// partition into parts of length sample_size
	for (unsigned short i = 0; i < all_paths.size(); i += sample_size) {
		auto last = min(all_paths.size(), (size_t) i+sample_size);
		vector<unsigned short> subset(all_paths.begin()+i, all_paths.begin()+last);
		sort(subset.begin(), subset.end());
		result.push_back(subset);
	}

	// if the last sample is shorter than the sample_size, randomly add additional paths
	unsigned short missing = sample_size - result.at(result.size()-1).size();
	if (missing > 0) {
		this->select_single_subset(result[result.size()-1], missing);
	}
}

void PathSampler::partition_samples(vector<vector<unsigned short>>& result, unsigned short sample_size) const {
	vector<pair<unsigned short,unsigned short>> all_samples;
	assert(this->total_number > 0);
	unsigned short n = this->total_number - 1;
	bool reference_added = this->total_number % 2 != 0;

	if (reference_added) {
		// reference not part of the panel
		for (size_t i = 1; i < n; i+=2) {
			all_samples.push_back(pair<size_t,size_t>(i,i+1));
		}
	} else {
		for (size_t i = 0; i < n; i+=2) {
			all_samples.push_back(pair<size_t,size_t>(i,i+1));
		}
	}

	random_shuffle(all_samples.begin(), all_samples.end());
	// list all paths
	vector<unsigned short> all_paths;

	if (reference_added) all_paths.push_back(0);

	for (auto p : all_samples) {
		all_paths.push_back(p.first);
		all_paths.push_back(p.second);
	}


	// partition into parts of length sample_size
	for (unsigned short i = 0; i < all_paths.size(); i += sample_size) {
		auto last = min(all_paths.size(), (size_t) i+sample_size);
		vector<unsigned short> subset(all_paths.begin()+i, all_paths.begin()+last);
		sort(subset.begin(), subset.end());
		result.push_back(subset);
	}
	// if the last sample is shorter than the sample_size, randomly add additional paths
	unsigned short missing;
	missing = sample_size - result.at(result.size()-1).size();

	if (missing > 0) {
		this->select_single_subset(result[result.size()-1], missing);
	}
}
