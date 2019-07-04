#include "catch.hpp"
#include "utils.hpp"
#include "../src/dnasequence.hpp"
#include <string>

using namespace std;

TEST_CASE("DnaSequence length", "[DnaSequence length]") {
	DnaSequence d;
	REQUIRE(d.length() == 0);
	string seq = "AAAATTGTG";
	d.append(seq);
	REQUIRE(d.length() == 9);

	string seq2 = "AAGTGTG";
	DnaSequence d2(seq2);
	REQUIRE(d2.length() == 7);
	string seq3 = "T";
	d2.append(seq3);
	REQUIRE(d2.length() == 8);
}

TEST_CASE("DnaSequence operator", "[DnaSequence operator]") {
	DnaSequence d;
	string seq = "ATGCTGGCTA";
	d.append(seq);
	for (size_t i = 0; i < d.length(); ++i) {
		REQUIRE(d[i] == seq[i]);
	}

	string seq2 = "NNNAAAA";
	d.append(seq2);
	REQUIRE(d.length() == 17);
	for (size_t i = 0; i < d.length(); ++i) {
		REQUIRE(d[i] == (seq+seq2)[i]);
	}
}

TEST_CASE("DnaSequence substr", "[DnaSequence substring]") {
	string seq = "ATGTTGTGCTGATGCTTGANNNTGATTCG";
	DnaSequence d(seq);
	REQUIRE(d.length() == seq.size());

	string substring;
	d.substr(3,17,substring);
	REQUIRE(substring == "TTGTGCTGATGCTT");

	substring.clear();
	d.substr(6,10,substring);
	REQUIRE(substring == "TGCT");

	substring.clear();
	d.substr(0,1,substring);
	REQUIRE(substring == "A");

	substring.clear();
	d.substr(0, d.length(),substring);
	REQUIRE(substring == seq);

	substring.clear();
	d.substr(0,0,substring);
	REQUIRE(substring == "");
}
