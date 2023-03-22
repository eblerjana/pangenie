#include "catch.hpp"
#include "utils.hpp"
#include "../src/kmerparser.hpp"
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("KmerParser empty_kmers", "[KmerParser empty_kmers]") {

	string line = "chr1\t1\2\t\t";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

}
