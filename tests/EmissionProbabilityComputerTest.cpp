#include "catch.hpp"
#include "utils.hpp"
#include "../src/emissionprobabilitycomputer.hpp"
#include "../src/copynumber.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("EmissionProbabilityComputer get_emission_probability", "EmissionProbabilityComputer [get_emission_probability]"){
	// construct UniqueKmers object
	vector<string> kmers = {"CATG", "ATGC", "CATT", "ATTG", "TTGC"};
	vector<vector<size_t>> paths = { {0}, {0}, {1,2}, {1,2}, {1,2} };
	vector<CopyNumber> cns = { CopyNumber(0.01, 0.2, 0.0), CopyNumber(0.001,0.5,0.001), CopyNumber(0.0,0.3,0.02), CopyNumber(0.05,0.6,0.0), CopyNumber(0.01,0.2,0.01)};
	UniqueKmers unique_kmers(0, 1000);
	for (unsigned int i = 0; i < kmers.size(); ++i) {
		unique_kmers.insert_kmer(cns[i],  paths[i]);
	}

	// construct EmissionProbabilityComputer
	EmissionProbabilityComputer emission_prob_comp (&unique_kmers);
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,0), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,1), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,0), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,2), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,0), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,2), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,1), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,1), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,2), 0.0));
}
