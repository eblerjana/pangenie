#include "catch.hpp"
#include "../src/pathsampler.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>

using namespace std;

TEST_CASE("PathSampler select_single_subset", "[PathSampler select_single_subset]") {
	PathSampler sampler (10);
	vector<unsigned short> sample = {};

	for (unsigned short i = 0; i < 30; ++i) {
		sample.clear();
		sampler.select_single_subset(sample, 7);
		REQUIRE(sample.size() == 7);

		// make sure each path was only selected once
		set<unsigned short> s(sample.begin(), sample.end());
		vector<unsigned short> singleton (s.begin(), s.end());
		sort(sample.begin(), sample.end());
		sort(singleton.begin(), singleton.end());
		REQUIRE(singleton == sample);

		// make sure all selected paths are contained in original set
		for (auto s : sample) {
			REQUIRE(s <= 10);
		}
	}
}


TEST_CASE("PathSampler select_multiple_subsets", "[PathSampler select_multiple_subsets]") {
	PathSampler sampler (20);
	vector<vector<unsigned short>> samples = {};

	for (unsigned short i = 0; i < 30; ++i) {
		samples.clear();
		sampler.select_multiple_subsets(samples, 7, 5);
		REQUIRE(samples.size() == 5);

		for (auto sample : samples) {
			// make sure each path was only selected once
			set<unsigned short> s(sample.begin(), sample.end());
			vector<unsigned short> singleton (s.begin(), s.end());
			sort(sample.begin(), sample.end());
			sort(singleton.begin(), singleton.end());
			REQUIRE(singleton == sample);
			// make sure all selected paths are contained in original set
			for (auto s : sample) {
				REQUIRE(s <= 20);
			}
		}
	}
}


TEST_CASE("PathSampler partition_paths", "[PathSampler partition_paths]") {
	PathSampler sampler (20);
	vector<vector<unsigned short>> partitions = {};

	sampler.partition_paths(partitions, 4);
	REQUIRE(partitions.size() == 5);
	
	vector<unsigned short> all_paths = {};
	for (auto v : partitions) {
		for (auto e : v) {
			all_paths.push_back(e);
		}
	}
	sort(all_paths.begin(), all_paths.end());

	vector<unsigned short> expected_paths = {};
	for (unsigned short i = 0; i < 20; ++i) {
		expected_paths.push_back(i);
	}

	REQUIRE(all_paths == expected_paths);
}

TEST_CASE("PathSampler partition_paths2", "[PathSampler partition_paths2]") {
	PathSampler sampler (25);
	vector<vector<unsigned short>> partitions = {};

	sampler.partition_paths(partitions, 4);
	REQUIRE(partitions.size() == 7);
	
	vector<unsigned short> all_paths = {};
	for (auto v : partitions) {
		for (auto e : v) {
			all_paths.push_back(e);
		}
	}
	sort(all_paths.begin(), all_paths.end());

	vector<unsigned short> expected_paths = {};
	for (unsigned short i = 0; i < 25; ++i) {
		expected_paths.push_back(i);
	}

	REQUIRE(includes(all_paths.begin(), all_paths.end(),
                  expected_paths.begin(), expected_paths.end()));
}

TEST_CASE("PathSampler partition_samples", "[PathSampler partition_samples]") {
	PathSampler sampler (20);
	vector<vector<unsigned short>> partitions = {};

	sampler.partition_samples(partitions, 4);
	REQUIRE(partitions.size() == 5);
	
	vector<unsigned short> all_paths = {};
	for (auto v : partitions) {
		for (auto e : v) {
			all_paths.push_back(e);
		}
	}
	sort(all_paths.begin(), all_paths.end());

	vector<unsigned short> expected_paths = {};
	for (unsigned short i = 0; i < 20; ++i) {
		expected_paths.push_back(i);
	}
	REQUIRE(all_paths == expected_paths);
}


TEST_CASE("PathSampler partition_samples2", "[PathSampler partition_samples]") {
	PathSampler sampler (9);
	vector<vector<unsigned short>> partitions = {};

	sampler.partition_samples(partitions, 5);
//	REQUIRE(partitions.size() == 2);
	
	vector<unsigned short> all_paths = {};
	for (auto v : partitions) {
		cout << " ----" << endl; 
		for (auto e : v) {
			all_paths.push_back(e);
			cout << e << endl;
		}
	}
	sort(all_paths.begin(), all_paths.end());

	vector<unsigned short> expected_paths = {0};
	for (unsigned short i = 0; i < 9; ++i) {
		expected_paths.push_back(i);
	}
	REQUIRE(all_paths == expected_paths);
}

