#include "catch.hpp"
#include "../src/sampledpanel.hpp"
#include <vector>
#include <string>
#include "utils.hpp"

using namespace std;

TEST_CASE("SampledPanel testcase1", "[SampledPanel testcase1]") {
	vector<int> path_to_allele = {0,1,2,0,3,2,1,1,0,0,0};
	SampledPanel s (path_to_allele, 0);
	for (size_t i = 0; i < path_to_allele.size(); ++i) {
		REQUIRE(s.get_allele_on_path(i) == path_to_allele[i]);
	}

	REQUIRE(path_to_allele == s.get_all_paths());
	REQUIRE(s.get_nr_paths() == path_to_allele.size());
	REQUIRE(s.get_unique_kmers() == 0);
}

TEST_CASE("SampledPanel get_specific_alleles", "[SampledPanel get_specific_alleles]") {

	vector<unsigned short> path_to_allele = {0,1,1,0,3,2,1,0,1,4};
	SampledPanel s (path_to_allele, 9);
	vector<unsigned short> specific_alleles = {1,2};
	vector<int> expected = {-1,0,0,-1,-1,1,0,-1,0,-1};

	REQUIRE(expected == s.get_specific_alleles(specific_alleles).get_all_paths());
	REQUIRE(s.get_unique_kmers() == 9);
}


TEST_CASE("SampledPanel get_specific_alleles2", "[SampledPanel get_specific_alleles2]") {

	vector<unsigned short> path_to_allele = {0,1,1,0,3,2,1,0,1,4};
	SampledPanel s (path_to_allele, 1);
	vector<unsigned short> specific_alleles = {};
	vector<int> expected = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

	REQUIRE(expected == s.get_specific_alleles(specific_alleles).get_all_paths());
	REQUIRE(s.get_unique_kmers() == 1);
}