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

TEST_CASE("FastaReader extract_name", "[FastaReader extract_name]") {
	FastaReader f("../tests/data/simple-fasta.fa");
	REQUIRE(f.get_size_of("chr01") == 1688);
	REQUIRE(f.get_size_of("chr02") == 2135);

	// try extracting non-existent chromosome
	REQUIRE_THROWS(f.extract_name("chrNone"));

	// extract existing chromosome
	vector<string> sequence_names_before;
	f.get_sequence_names(sequence_names_before);

	vector<string> expected = {"chr01", "chr02"};
	REQUIRE(sequence_names_before == expected);
	FastaReader extracted = f.extract_name("chr01");

	// chr02 should still be there, but not chr01
	REQUIRE(f.get_size_of("chr02") == 2135);
	REQUIRE_THROWS(f.get_size_of("chr01"));

	// chr01 should now only be present in new FastaReader
	REQUIRE(extracted.contains_name("chr01"));
	REQUIRE(!f.contains_name("chr01"));

	vector<string> sequence_names_after;
	f.get_sequence_names(sequence_names_after);
	expected = {"chr02"};
	REQUIRE(sequence_names_after == expected);
	sequence_names_after.clear();
	extracted.get_sequence_names(sequence_names_after);
	expected = {"chr01"};
	REQUIRE(sequence_names_after == expected);

	REQUIRE(extracted.get_size_of("chr01") == 1688);

	// subsequence can no longer be extracted from non-existent variant
	string sequence;
	REQUIRE_THROWS(f.get_subsequence("chr01", 0, 10, sequence));

	// same sequence cannot be extracted twice
	REQUIRE_THROWS(f.extract_name("chr01"));
}

TEST_CASE("FastaReader extract_name2", "[FastaReader extract_name2]") {
	FastaReader f("../tests/data/simple-fasta.fa");
	FastaReader extracted;

	f.extract_name("chr01");
	f.extract_name("chr02");

	// remove all sequences
	REQUIRE(!f.contains_name("chr01"));
	REQUIRE(!f.contains_name("chr02"));

	// "f" should now no longer contain any sequences
	vector<string> names;
	f.get_sequence_names(names);
	REQUIRE(names.empty());
}