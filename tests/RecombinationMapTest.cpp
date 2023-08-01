#include "catch.hpp"
#include "utils.hpp"
#include "../src/recombinationmap.hpp"
#include "../src/variantreader.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("RecombinationMap load_genetic_map", "[RecombinationMap load_genetic_map]"){
	string filename = "../tests/data/recomb_map.map";
	vector<pair<size_t, long double>> g;
	load_genetic_map(filename, &g);
	REQUIRE(g.size() == 7);
}

TEST_CASE("RecombinationMap parse_map_inputs", "[RecombinationMap parse_map_inputs]") {
	string filename = "../tests/data/recomb_map_config.tsv";
	map<string, string> result;
	parse_map_inputs(filename, result);

	REQUIRE(result.size() == 4);
	REQUIRE(result["chr1"] == "filename1");
	REQUIRE(result["chr2"] == "filename2");
	REQUIRE(result["chr3"] == "filename3");
	REQUIRE(result["chr4"] == "filename4");
}

TEST_CASE("RecombinationMap load_genetic_map_malformed", "[RecombinationMap load_genetic_map_malformed]"){
	string filename = "../tests/data/recomb_map_malformed.map";
	vector<pair<size_t, long double>> g;
	REQUIRE_THROWS(load_genetic_map(filename, &g));
}

TEST_CASE("RecombinationMap compute_recombination_probability_uniform", "[RecombinationMap compute_recombination_probability_uniform]") {
	RecombinationMap r (1.26);
	REQUIRE(doubles_equal(r.compute_recombination_probability(1000, 1, 3500, 2), 0.00315));
	REQUIRE_THROWS(r.compute_recombination_probability(4000, 1, 3500, 2));
	REQUIRE_THROWS(r.compute_recombination_probability(1000, 3, 3500, 1));
}

TEST_CASE("RecombinationMap compute_recombination_probability_map", "RecombinationMap compute_recombination_probability_map") {
	string vcf = "../tests/data/small1.vcf";
	string fasta = "../tests/data/small1.fa";
	string filename = "../tests/data/recomb_map2.map";
	VariantReader v(vcf, fasta, 10, true);
	RecombinationMap r(filename, &v, "chrA");

	vector<double> expected = {2.8442017983318973e-05, 5.688403596663811e-06, 0.0001934057222865691, 5.688403596663792e-05, 5.7452876326304374e-05, 5.745287632630432e-05};
	vector<double> computed = {};

	vector<size_t> positions = {101, 151, 161, 501, 606, 702, 803};
	for (size_t i = 0; i < positions.size()-1; ++i) {
		computed.push_back(r.compute_recombination_probability(positions[i], i, positions[i+1], i+1));
	}

	REQUIRE(compare_vectors(computed, expected));
}

TEST_CASE("RecombinationMap compute_recombination_probability_map2", "RecombinationMap compute_recombination_probability_map2") {
	string vcf = "../tests/data/small3.vcf";
	string fasta = "../tests/data/small3.fa";
	string filename = "../tests/data/recomb_map3.map";
	VariantReader v(vcf, fasta, 10, true);
	RecombinationMap r(filename, &v, "chr1");

	vector<double> expected = {0.0001328663819544218, 7.247257197517953e-05, 0.00022043740642430087, 3.925597648635204e-05, 9.059071496908544e-05, 3.925597648635204e-05, 0.00031404781189237063, 0.00010568916746356294, 9.663009596683203e-05, 0.0005888396472979451, 8.757102447010112e-05, 0.00024459493041595337, 9.361040546784771e-05};
	vector<double> computed = {};

	vector<double> positions = {};

	for (size_t i = 0; i < v.size_of("chr1"); i++) {
		positions.push_back(v.get_variant("chr1", i).get_start_position());
	}

	for (size_t i = 0; i < positions.size()-1; ++i) {
		computed.push_back(r.compute_recombination_probability(positions[i], i, positions[i+1], i+1));
	}

	REQUIRE(compare_vectors(computed, expected));
}

TEST_CASE("RecombinationMap compute_recombination_probability_map3", "RecombinationMap compute_recombination_probability_map3") {
	string vcf = "../tests/data/small3.vcf";
	string fasta = "../tests/data/small3.fa";
	string filename = "../tests/data/recomb_map4.map";
	VariantReader v(vcf, fasta, 10, true);
	RecombinationMap r(filename, &v, "chr1");

	vector<double> expected = {6.778756507563877e-05, 3.697503549571213e-05, 0.00011246573296641582, 2.0028144226857947e-05, 4.6218794369723426e-05, 2.002814422696897e-05, 0.00016022515381519664, 5.392192676456631e-05, 4.9300047327838215e-05, 0.0003004221634033133, 4.4678167890888076e-05, 0.00012479074479831986, 4.77594208486698e-05};
	vector<double> computed = {};

	vector<double> positions = {};

	for (size_t i = 0; i < v.size_of("chr1"); i++) {
		positions.push_back(v.get_variant("chr1", i).get_start_position());
	}

	for (size_t i = 0; i < positions.size()-1; ++i) {
		computed.push_back(r.compute_recombination_probability(positions[i], i, positions[i+1], i+1));
	}

	REQUIRE(compare_vectors(computed, expected));
}