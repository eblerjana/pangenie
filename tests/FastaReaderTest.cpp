#include "catch.hpp"
#include "utils.hpp"
#include "../src/fastareader.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("FastaReader contains_name", "[FastaReader contains_name]") {
	FastaReader f("../tests/data/simple-fasta.fa");
	REQUIRE(f.contains_name("chr01"));
	REQUIRE(f.contains_name("chr02"));
	REQUIRE(!f.contains_name("chr03"));
}

TEST_CASE("FastaReader get_size_of", "[FastaReader get_size_of]") {
	FastaReader f("../tests/data/simple-fasta.fa");
	REQUIRE(f.get_size_of("chr01") == 1688);
	REQUIRE(f.get_size_of("chr02") == 2135);
	REQUIRE_THROWS(f.get_size_of("chrNone"));
	REQUIRE(f.get_total_kmers(20) == 3785);
}

TEST_CASE("FastaReader get_subsequence", "[FastaReader get_subsequence]") {
	FastaReader f("../tests/data/simple-fasta.fa");
	string sequence;
	f.get_subsequence("chr01", 0, 10, sequence);
	REQUIRE(sequence == "CATTTTAAAG");

	sequence.clear();
	f.get_subsequence("chr01", 21, 40, sequence);
	REQUIRE(sequence == "CCCAGAGCAGGCAAAACCC");

	sequence.clear();
	f.get_subsequence("chr02", 1, 12, sequence);
	REQUIRE(sequence == "CCAACAATTTA");

	sequence.clear();
	f.get_subsequence("chr02", 71, 81, sequence);
	REQUIRE(sequence == "TCAAATCACA");

	sequence.clear();
	REQUIRE_THROWS(f.get_subsequence("chrNone", 71, 80, sequence));
}

TEST_CASE("FastaReader invalid", "[FastaReader invalid]") {
	REQUIRE_THROWS(FastaReader("../tests/data/broken-fasta.fa"));
}
