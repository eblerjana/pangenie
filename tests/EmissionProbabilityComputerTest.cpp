#include "catch.hpp"
#include "../src/emissionprobabilitycomputer.h"
#include "../src/copynumber.h"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("EmissionProbabilityComputer simple example", "EmissionProbabilityComputer [simple]"){
	vector<string> kmers = {"ATGCA", "TTGAC", "GGTGA"};
	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1),CopyNumber(0.1,0.2,0.7)};
	vector<vector<int>> paths = { {0,1}, {0,2}, {2,3} };
	EmissionProbabilityComputer e(4);
	for (size_t i = 0; i < 3; ++i){
		e.insert_kmer(kmers[i], cns[i], paths[i]);
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(e.contains_kmer(kmers[i]));
		REQUIRE(e.get_copynumber_of(kmers[i]) == cns[i]);
		REQUIRE(e.get_paths_of(kmers[i]) == paths[i]);
	}
}

TEST_CASE("EmissionProbabilityComputer check exception", "[EmissionProbabilityComputer exception]"){
	EmissionProbabilityComputer e(3);
	vector<string> kmers = {"ATGCA", "TTGAC", "GGTGA"};
	for (size_t i = 0; i < 3; ++i){
		CHECK_THROWS(e.get_copynumber_of(kmers[i]));
		CHECK_THROWS(e.get_paths_of(kmers[i]));
	}
	CHECK_THROWS(e.insert_kmer("ATTGA", CopyNumber(0.1,0.5,0.4), {3}));
}

TEST_CASE("EmissionProbabiliyComputer emission", "[EmissionProbabilityComputer emission]"){
	EmissionProbabilityComputer e(3);
	vector<string> kmers = {"ATGCA", "TTGAC", "GGTGA"};
	vector<CopyNumber> cns = {CopyNumber(0.1,0.8,0.1), CopyNumber(0.2,0.7,0.1), CopyNumber(0.2,0.7,0.1)};
	vector<vector<int>> paths = { {0}, {0}, {1,2} };
	for (size_t i = 0; i < 3; ++i){
		e.insert_kmer(kmers[i], cns[i], paths[i]);
	}
	// long double get_emission_probability(int path1, int path2);
	
}
