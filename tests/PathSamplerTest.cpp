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
	vector<size_t> sample = {};

	for (size_t i = 0; i < 30; ++i) {
		sample.clear();
		sampler.select_single_subset(sample, 7);
		REQUIRE(sample.size() == 7);

		// make sure each path was only selected once
		set<size_t> s(sample.begin(), sample.end());
		vector<size_t> singleton (s.begin(), s.end());
		sort(sample.begin(), sample.end());
		sort(singleton.begin(), singleton.end());
		REQUIRE(singleton == sample);

		// make sure all selected paths are contained in original set
		for (auto s : sample) {
			REQUIRE(s <= 10);
			REQUIRE(s >= 0);
		}
	}
}


TEST_CASE("PathSampler select_multiple_subsets", "[PathSampler select_multiple_subsets]") {
	PathSampler sampler (20);
	vector<vector<size_t>> samples = {};

	for (size_t i = 0; i < 30; ++i) {
		samples.clear();
		sampler.select_multiple_subsets(samples, 7, 5);
		REQUIRE(samples.size() == 5);

		for (auto sample : samples) {
			// make sure each path was only selected once
			set<size_t> s(sample.begin(), sample.end());
			vector<size_t> singleton (s.begin(), s.end());
			sort(sample.begin(), sample.end());
			sort(singleton.begin(), singleton.end());
			REQUIRE(singleton == sample);
			// make sure all selected paths are contained in original set
			for (auto s : sample) {
				REQUIRE(s <= 20);
				REQUIRE(s >= 0);
			}
		}
	}
}
