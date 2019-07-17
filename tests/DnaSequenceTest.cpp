#include "catch.hpp"
#include "utils.hpp"
#include "../src/dnasequence.hpp"
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("DnaSequence size", "[DnaSequence length]") {
	DnaSequence d;
	REQUIRE(d.size() == 0);
	string seq = "AAAATTGTG";
	d.append(seq);
	REQUIRE(d.size() == 9);

	string seq2 = "AAGTGTG";
	DnaSequence d2(seq2);
	REQUIRE(d2.size() == 7);
	string seq3 = "T";
	d2.append(seq3);
	REQUIRE(d2.size() == 8);
}

TEST_CASE("DnaSequence operator", "[DnaSequence operator]") {
	DnaSequence d;
	string seq = "ATGCTGGCTA";
	d.append(seq);
	for (size_t i = 0; i < d.size(); ++i) {
		REQUIRE(d[i] == seq[i]);
	}

	string seq2 = "NNNAAAA";
	d.append(seq2);
	REQUIRE(d.size() == 17);
	for (size_t i = 0; i < d.size(); ++i) {
		REQUIRE(d[i] == (seq+seq2)[i]);
	}
}

TEST_CASE("DnaSequence substr (string)", "[DnaSequence substr (string)]") {
	string seq = "ATGTTGTGCTGATGCTTGANNNTGATTCG";
	DnaSequence d(seq);
	REQUIRE(d.size() == seq.size());

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
	d.substr(0, d.size(),substring);
	REQUIRE(substring == seq);

	substring.clear();
	d.substr(0,0,substring);
	REQUIRE(substring == "");
}

TEST_CASE("DnaSequence substr (DnaSequence)", "[DnaSequence substr (DnaSequence)]") {
	string seq = "ATGTTGTGCTGATGCTTGANNNTGATTCG";
	DnaSequence d(seq);

	DnaSequence substring;
	d.substr(3, 17, substring);
	REQUIRE(substring.to_string() == "TTGTGCTGATGCTT");
	REQUIRE(substring.size() == 14);

	substring.clear();
	d.substr(6, 10, substring);
	REQUIRE(substring.to_string() == "TGCT");
	REQUIRE(substring.size() == 4);

	substring.clear();
	d.substr(0, 1, substring);
	REQUIRE(substring.to_string() == "A");
	REQUIRE(substring.size() == 1);

	substring.clear();
	d.substr(0, d.size(), substring);
	REQUIRE(substring.to_string() == seq);
	REQUIRE(substring.size() == seq.size());

	substring.clear();
	d.substr(0, 0, substring);
	REQUIRE(substring.to_string() == "");
	REQUIRE(substring.size() == 0);
}

TEST_CASE("DnaSequence reverse", "[DnaSequence reverse]") {
	string seq = "AGGTGCCTGATCGTACGGCTAGCTGAACNNNT";
	string rev_seq = "TNNNCAAGTCGATCGGCATGCTAGTCCGTGGA";
	DnaSequence d(seq);
	d.reverse();
	REQUIRE(d.to_string() == rev_seq);
}

TEST_CASE("DnaSequence reverse2", "[DnaSequence reverse2]") {
	string seq = "ATTGC";
	string rev_seq = "CGTTA";
	DnaSequence d(seq);
	d.reverse();
	REQUIRE(d.to_string() == rev_seq);
}

TEST_CASE("DnaSequence reverse_complement", "DnaSequence reverse_complement") {
	string seq = "AGTGAAACGCGTTCGAATGCNNNTGCGNA";
	string rev_compl_seq = "TNCGCANNNGCATTCGAACGCGTTTCACT";
	DnaSequence d(seq);
	d.reverse_complement();
	REQUIRE(d.to_string() == rev_compl_seq);
}

TEST_CASE("DnaSequence reverse_complement2", "[DnaSequence reverse_complement2]") {
	string seq = "ATGCG";
	string rev_compl_seq = "CGCAT";
	DnaSequence d(seq);
	d.reverse_complement();
	REQUIRE(d.to_string() == rev_compl_seq);
}

TEST_CASE("DnaSequence append (DnaSequence)", "[DnaSequence append (DnaSequence)]") {
	string seq = "AAAAAAA";
	string seq2 = "TTTTT";
	string expected = seq + seq2;
	DnaSequence d(seq);
	DnaSequence d2(seq2);
	d.append(d2);
	REQUIRE(d.to_string() == expected);

	vector<string> seqs = {"AAA", "ATGC", "T", "A", "CCC"};
	vector<DnaSequence> dna_seqs;
	for (auto s : seqs) {
		dna_seqs.push_back(DnaSequence(s));
	}

	DnaSequence result;
	for (auto d : dna_seqs) {
		result.append(d);
	}
	REQUIRE(result.size() == 12);
	REQUIRE(result.to_string() == "AAAATGCTACCC");
}

TEST_CASE("DnaSequence operators", "[DnaSequence operators]") {
	string seq = "ATGGTGT";
	string seq2 = "ATGTGT";
	string seq3 = "ATGGTGT";
	DnaSequence d(seq);
	DnaSequence d2(seq2);
	REQUIRE(d != d2);
	DnaSequence d3(seq3);
	REQUIRE(d == d3);
}

TEST_CASE("DnaSequence base_at", "[DnaSequence base_at]") {
	string seq = "ATGGTC";
	DnaSequence d(seq);
	DnaSequence base = d.base_at(1);
	REQUIRE(base.to_string() == "T");
	REQUIRE(base.size() == 1);
	base = d.base_at(4);
	REQUIRE(base.to_string() == "T");
	REQUIRE(base.size() == 1);
}
