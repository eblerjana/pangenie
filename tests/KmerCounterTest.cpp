#include "catch.hpp"
#include "utils.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/jellyfishreader.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("JellyfishCounter", "[JellyfishCounter]") {
	JellyfishCounter counter("../tests/data/reads.fa", 10);
	string read = "ATGCTGTAAAAAAACGGC";
	for (size_t i = 0; i < read.size()-9; ++i) {
		string kmer = read.substr(i,10);
		REQUIRE(counter.getKmerAbundance(kmer) == 1);
	}
}

TEST_CASE("JellyfishReader", "[JellyfishReader]") {
	JellyfishReader reader ("../tests/data/reads.jf", 10);
	string read = "ATGCTGTAAAAAAACGGC";
	for (size_t i = 0; i < read.size()-9; ++i) {
		string kmer = read.substr(i,10);
		REQUIRE(reader.getKmerAbundance(kmer) == 1);
	}

	// reads where counted without the -C option
	REQUIRE_THROWS(JellyfishReader("../tests/data/reads.no-canonical.jf", 10));

	// wrong kmer size used
	REQUIRE_THROWS(JellyfishReader("../tests/data/reads.jf", 11));

}
