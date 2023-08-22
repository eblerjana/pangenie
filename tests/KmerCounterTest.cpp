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
	KmerCounter* t = new JellyfishCounter ("../tests/data/reads.fa", 10);
	delete t;
}

TEST_CASE("JellyfishCounter_if", "[JellyfishCounter_if]") {
	JellyfishCounter counter("../tests/data/reads.fa", {"../tests/data/kmerfile.fa"}, 10);
	// these two kmers are in kmerfile.fa and should have been counted
	REQUIRE(counter.getKmerAbundance("ATGCTGTAAA") == 1);
	REQUIRE(counter.getKmerAbundance("TGCTGTAAAA") == 1);
	// the following kmers are not contained in kmerfile.fa, thus they should have count 0
	string kmers = "GCTGTAAAAAAACGGC";
	for (size_t i = 0; i < kmers.size()-9; ++i) {
		string kmer = kmers.substr(i,10);
		REQUIRE(counter.getKmerAbundance(kmer) == 0);
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
