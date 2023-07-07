#include "catch.hpp"
#include "utils.hpp"
#include "../src/recombinationmap.hpp"
#include "../src/variantreader.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("RecombinationMap load_genetic_map", "[RecombinationMap load_genetic_map]"){
	string filename = "../tests/data/recomb_map.map";
	vector<pair<size_t, float>> g;
	load_genetic_map(filename, &g);
	REQUIRE(g.size() == 7);
}

TEST_CASE("RecombinationMap load_genetic_map_malformed", "[RecombinationMap load_genetic_map_malformed]"){
	string filename = "../tests/data/recomb_map_malformed.map";
	vector<pair<size_t, float>> g;
	REQUIRE_THROWS(load_genetic_map(filename, &g));
}

TEST_CASE("RecombinationMap compute_recombination_probability_uniform", "[RecombinationMap compute_recombination_probability_uniform]") {
	RecombinationMap r (1.26);
	REQUIRE(doubles_equal(r.compute_recombination_probability(1000, 1, 3500, 2), 0.00315));
	REQUIRE_THROWS(doubles_equal(r.compute_recombination_probability(4000, 1, 3500, 2), 0.00315));
	REQUIRE_THROWS(doubles_equal(r.compute_recombination_probability(1000, 1, 3500, 3), 0.00315));
}

TEST_CASE("RecombinationMap compute_recombination_probability_map", "RecombinationMap compute_recombination_probability_map") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	string filename = "../tests/data/recomb_map.map";
	VariantReader v(vcf, fasta, 10, true);
	RecombinationMap r(filename, &v, "chrA");

	vector<size_t> positions = {101, 151, 161, 166, 501, 606, 702};
	for (size_t i = 0; i < positions.size()-1; ++i) {
		cout << r.compute_recombination_probability(positions[i], i, positions[i+1], i+1) << endl;
	}
}