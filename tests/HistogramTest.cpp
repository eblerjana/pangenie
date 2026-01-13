#include "catch.hpp"
#include "utils.hpp"
#include "../src/sequenceutils.hpp"
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

	for (auto p : peak_values) {
		cout << p << " ";
	}
	cout << endl;

	REQUIRE(peak_ids.size() == 1);
	REQUIRE(peak_values.size() == 1);
	REQUIRE(peak_ids[0] == 1);
	REQUIRE(peak_values[0] == 4);
}

TEST_CASE("Histogram2", "[Histogram2]") {
	Histogram histo("../tests/data/test.histo", 10000);
	histo.smooth_histogram();
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histo.find_peaks(peak_ids, peak_values);
	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, true);
	REQUIRE(kmer_coverage_estimate == 56);
}

TEST_CASE("Histogram3", "[Histogram3]") {
	Histogram histo("../tests/data/test2.histo", 10000);
	histo.smooth_histogram();
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histo.find_peaks(peak_ids, peak_values);
	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, true);
	REQUIRE(kmer_coverage_estimate == 26);
}


TEST_CASE("Histogram4", "[Histogram4]") {
	Histogram histo("../tests/data/test3.histo", 10000);
	histo.smooth_histogram();
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histo.find_peaks(peak_ids, peak_values);
	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, true);
	REQUIRE(kmer_coverage_estimate == 60);
}

TEST_CASE("Histogram5", "[Histogram5]") {
	Histogram histo("../tests/data/test4.histo", 10000);
	histo.smooth_histogram();
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histo.find_peaks(peak_ids, peak_values);
	size_t kmer_coverage_estimate = compute_kmer_coverage(peak_ids, peak_values, true);
	REQUIRE(kmer_coverage_estimate == 42);
}
