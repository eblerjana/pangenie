#include "catch.hpp"
#include "utils.hpp"
#include "../src/emissionprobabilitycomputer.hpp"
#include "../src/copynumber.hpp"
#include "../src/probabilitytable.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace std;

TEST_CASE("EmissionProbabilityComputer get_emission_probability", "EmissionProbabilityComputer [get_emission_probability]"){
	// construct UniqueKmers object
	vector<string> kmers = {"CATG", "ATGC", "CATT", "ATTG", "TTGC"};
	vector<vector<unsigned char>> alleles = {{0}, {0}, {1}, {1}, {1}};
	vector<unsigned short> counts = {4, 6, 8, 2, 5};
	vector<CopyNumber> cns = { CopyNumber(0.01, 0.2, 0.0), CopyNumber(0.001,0.5,0.001), CopyNumber(0.0,0.3,0.02), CopyNumber(0.05,0.6,0.0), CopyNumber(0.01,0.2,0.01)};
	ProbabilityTable probs (0,10,10,0.0);
	vector<unsigned char> path_to_allele = {0, 1, 1};
	shared_ptr<UniqueKmers> unique_kmers = shared_ptr<UniqueKmers> (new UniqueKmers(1000, path_to_allele));
	for (unsigned int i = 0; i < kmers.size(); ++i) {
		unique_kmers->insert_kmer(counts[i],  alleles[i]);
		probs.modify_probability(0, counts[i], cns[i]);
	}
	vector<unsigned short> path_ids;
	vector<unsigned char> allele_ids;
	unique_kmers->get_path_ids(path_ids, allele_ids);
	// construct EmissionProbabilityComputer
	EmissionProbabilityComputer emission_prob_comp (unique_kmers, &probs);
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,0), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,1), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,0), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,1), 0.0));
}


TEST_CASE("EmissionProbabilityComputer get_emission_probability_undefined1", "EmissionProbabilityComputer [get_emission_probability_undefined1]"){
	// construct UniqueKmers object
	vector<string> kmers = {"CATG", "ATGC", "CATT", "ATTG", "TTGC"};
	vector<vector<unsigned char>> alleles = {{0}, {0}, {1}, {1}, {1}};
	vector<unsigned short> counts = {4, 6, 8, 2, 5};
	vector<CopyNumber> cns = { CopyNumber(0.01, 0.2, 0.0), CopyNumber(0.001,0.5,0.001), CopyNumber(0.0,0.3,0.02), CopyNumber(0.05,0.6,0.0), CopyNumber(0.01,0.2,0.01)};
	ProbabilityTable probs (0,10,10,0.0);
	vector<unsigned char> path_to_allele = {0, 1, 2};
	shared_ptr<UniqueKmers> unique_kmers = shared_ptr<UniqueKmers> (new UniqueKmers(1000, path_to_allele));
	unique_kmers->set_undefined_allele(2);
	for (unsigned int i = 0; i < kmers.size(); ++i) {
		unique_kmers->insert_kmer(counts[i],  alleles[i]);
		probs.modify_probability(0, counts[i], cns[i]);
	}

	vector<unsigned short> path_ids;
	vector<unsigned char> allele_ids;
	unique_kmers->get_path_ids(path_ids, allele_ids);

	// construct EmissionProbabilityComputer
	EmissionProbabilityComputer emission_prob_comp (unique_kmers, &probs);
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,0), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,1), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,0), 0.0036));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,1), 0.0));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(0,2), 0.000128225));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,0), 0.000128225));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(1,2), 0.000132565));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,1), 0.000132565));
	REQUIRE (doubles_equal(emission_prob_comp.get_emission_probability(2,2), 0.000019852));
	
}
