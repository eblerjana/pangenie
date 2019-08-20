#include "catch.hpp"
#include "utils.hpp"
#include "../src/kmercounter.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("KmerCounter", "[KmerCounter]") {
	KmerCounter counter("../tests/data/reads.fa", 10);
	string read = "ATGCTGTAAAAAAACGGC";
	for (size_t i = 0; i < read.size()-9; ++i) {
		string kmer = read.substr(i,10);
		REQUIRE(counter.getKmerAbundance(kmer) == 1);
	}
}
