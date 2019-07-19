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
	vector<vector<unsigned char>> alleles = {{0}, {0}, {1}, {1}, {1}};
	vector<CopyNumber> cns = { CopyNumber(0.01, 0.2, 0.0), CopyNumber(0.001,0.5,0.001), CopyNumber(0.0,0.3,0.02), CopyNumber(0.05,0.6,0.0), CopyNumber(0.01,0.2,0.01)};
	UniqueKmers unique_kmers(0, 1000);
	unique_kmers.insert_path(0,0);
	unique_kmers.insert_path(1,1);
	unique_kmers.insert_path(2,1);
	for (unsigned int i = 0; i < kmers.size(); ++i) {
		unique_kmers.insert_kmer(cns[i],  alleles[i]);
	}

	// construct ColumnIndexer
	ColumnIndexer column_indexer(0,9);

	vector<size_t> path_ids;
	vector<unsigned char> allele_ids;
	unique_kmers.get_path_ids(path_ids, allele_ids);

	// insert states
	for (size_t i = 0; i < path_ids.size(); ++i) {
		for (size_t j = 0; j < path_ids.size(); ++j) {
			column_indexer.insert( make_pair(path_ids[i],path_ids[j]), make_pair(allele_ids[i], allele_ids[j]));
		}
	}

	// construct EmissionProbabilityComputer
	EmissionProbabilityComputer emission_prob_comp (&unique_kmers, &column_indexer);
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(3), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(6), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(5), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(7), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(4), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(8), 0.0));
}
