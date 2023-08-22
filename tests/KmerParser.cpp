#include "catch.hpp"
#include "utils.hpp"
#include "../src/kmerparser.hpp"
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("KmerParser all_empty_kmers", "[KmerParser all_empty_kmers]") {

	string line = "chr1\t1\t2\tnan\tnan";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

	vector<string> kmers;
	vector<string> flanking_kmers;
	string chrom;
	size_t start;
	bool is_header;
	parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);
	REQUIRE(kmers.size() == 0);
	REQUIRE(flanking_kmers.size() == 0);
	REQUIRE(chrom == "chr1");
	REQUIRE(start == 1);
}


TEST_CASE("KmerParser empty_unique_kmers", "[KmerParser empty_unique_kmers]") {

	string line = "chr1\t1\t2\tnan\tAAAT,TGGG";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

	vector<string> kmers;
	vector<string> flanking_kmers;
	string chrom;
	size_t start;
	bool is_header;
	parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);
	REQUIRE(kmers.size() == 0);
	REQUIRE(flanking_kmers.size() == 2);
	REQUIRE(flanking_kmers[0] == "AAAT");
	REQUIRE(flanking_kmers[1] == "TGGG");
	REQUIRE(chrom == "chr1");
	REQUIRE(start == 1);
}


TEST_CASE("KmerParser empty_flanking_kmers", "[KmerParser empty_flanking_kmers]") {

	string line = "chr1\t1\t2\tTGTG,ATGT\tnan";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

	vector<string> kmers;
	vector<string> flanking_kmers;
	string chrom;
	size_t start;
	bool is_header;
	parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);
	REQUIRE(kmers.size() == 2);
	REQUIRE(kmers[0] == "TGTG");
	REQUIRE(kmers[1] == "ATGT");
	REQUIRE(flanking_kmers.size() == 0);
	REQUIRE(chrom == "chr1");
	REQUIRE(start == 1);
}


TEST_CASE("KmerParser parse_kmers", "[KmerParser parse_kmers]") {

	string line = "chr1\t1\t2\tTGTG,ATGT\tTTTT,GGGG";
	vector<string> result;
	parse(result, line, '\t');
	REQUIRE(result.size() == 5);

	vector<string> kmers;
	vector<string> flanking_kmers;
	string chrom;
	size_t start;
	bool is_header;
	parse_kmer_line(line, chrom, start, kmers, flanking_kmers, is_header);
	REQUIRE(kmers.size() == 2);
	REQUIRE(kmers[0] == "TGTG");
	REQUIRE(kmers[1] == "ATGT");
	REQUIRE(flanking_kmers.size() == 2);
	REQUIRE(flanking_kmers[0] == "TTTT");
	REQUIRE(flanking_kmers[1] == "GGGG");
	REQUIRE(chrom == "chr1");
	REQUIRE(start == 1);
}