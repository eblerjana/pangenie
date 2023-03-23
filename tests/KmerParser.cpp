#include "catch.hpp"
#include "utils.hpp"
#include "../src/kmerparser.hpp"
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("KmerParser empty_kmers", "[KmerParser empty_kmers]") {

	string line = "chr1\t1\t2\tnan\tnan";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

	vector<string> kmers;
	vector<string> flanking_kmers;
	bool is_header;
	parse_kmer_line(line, kmers, flanking_kmers, is_header);
	REQUIRE(kmers.size() == 0);
	REQUIRE(flanking_kmers.size() == 0);

}
