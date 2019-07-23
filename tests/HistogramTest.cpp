#include "catch.hpp"
#include "utils.hpp"
#include "../src/histogram.hpp"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("Histogram", "[Histogram]") {
	Histogram histo(10);
	vector<size_t> values = {0,0,1,1,1,1,2,2,3};
	for (auto v : values) {
		histo.add_value(v);
	}	
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histo.find_peaks(peak_ids, peak_values);
	REQUIRE(peak_ids.size() == 1);
	REQUIRE(peak_values.size() == 1);
	REQUIRE(peak_ids[0] == 1);
	REQUIRE(peak_values[0] == 4);
}
