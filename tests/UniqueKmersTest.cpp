#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("UniqueKmers testcase 1", "[UniqueKmers testcase 1]"){
	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1),CopyNumber(0.1,0.2,0.7)};
	vector<vector<size_t>> paths = {{0,1,2}, {0,1}, {2}};
	UniqueKmers u(0,1000);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(cns[i], alleles[i]);
		for (auto& p : paths[i]){
			REQUIRE(u.kmer_on_path(i, p));
		}
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(u.get_variant_index() == 0);
		REQUIRE(u.get_variant_position() == 1000);
		REQUIRE(u.get_copynumber_of(i) == cns[i]);
	}
}

TEST_CASE("UniqueKmers get_copynumber_of", "[UniqueKmers get_copynumber_of]"){
	UniqueKmers u(1, 2000);
	for (size_t i = 0; i < 3; ++i){
		CHECK_THROWS(u.get_copynumber_of(i));
	}
	CopyNumber cn(0.05, 0.9, 0.05);
	vector<unsigned char> alleles = {0};
	u.insert_kmer(cn, alleles);
	REQUIRE(u.get_copynumber_of(0) == cn);
}

TEST_CASE("UniqueKmers insert_empty_path", "[UniqueKmers insert_empty_path]") {
	UniqueKmers u(0, 1000);
	vector<CopyNumber> cns = { CopyNumber(0.001, 0.5, 0.001), CopyNumber(1.0, 0.0, 0.0), CopyNumber(0.001,0.001, 0.6) };
	vector<vector<size_t>> paths = { {0}, {0}, {2} };
	vector<vector<unsigned char>> alleles = {{0}, {0}, {1}};

	u.insert_path(0,0);
	u.insert_path(2,1);

	// insert the kmers
	for (size_t i = 0; i < 3; ++i) {
		u.insert_kmer(cns[i], alleles[i]);
		for (auto& p : paths[i]) {
			REQUIRE(u.kmer_on_path(i, p));
		}
	}

	vector<size_t> path_ids;
	vector<unsigned char> allele_ids;
	u.get_path_ids(path_ids, allele_ids);
	REQUIRE(u.size() == 3);

	REQUIRE(path_ids.size() == 2);
	REQUIRE(path_ids[0] == 0);
	REQUIRE(path_ids[1] == 2);

	REQUIRE(allele_ids.size() == 2);
	REQUIRE(allele_ids[0] == 0);
	REQUIRE(allele_ids[1] == 1);

	// insert an empty path
	u.insert_empty_allele(2);
	u.insert_path(1,2);
	path_ids.clear();
	allele_ids.clear();
	u.get_path_ids(path_ids, allele_ids);
	REQUIRE(u.size() == 3);

	REQUIRE(path_ids.size() == 3);
	REQUIRE(path_ids[0] == 0);
	REQUIRE(path_ids[1] == 1);
	REQUIRE(path_ids[2] == 2);

	REQUIRE(allele_ids.size() == 3);
	REQUIRE(allele_ids[0] == 0);
	REQUIRE(allele_ids[1] == 2);
	REQUIRE(allele_ids[2] == 1);

	// make sure path is indeed empty
	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(! u.kmer_on_path(i, 1));
	}

	// add kmer to the empty path
	vector<unsigned char> allele = {2};
	u.insert_kmer(CopyNumber(0.01, 0.8, 0.01), allele);
	REQUIRE(u.size() == 4);
	REQUIRE(u.kmer_on_path(3, 1));
}

TEST_CASE("UniqueKmers insert_empty_path2", "[UniqueKmers insert_empty_path2]") {
	UniqueKmers u(0, 1000);
	CopyNumber cn(0.9, 0.1, 0.2);
	vector<unsigned char> allele = {1};

	u.insert_kmer(cn, allele);
	u.insert_path(1,1);
	REQUIRE(u.size() == 1);
	REQUIRE(u.kmer_on_path(0, 1));

	vector<size_t> path_ids;
	vector<unsigned char> allele_ids;
	u.get_path_ids(path_ids, allele_ids);
	REQUIRE(path_ids.size() == 1);
	REQUIRE(path_ids[0] == 1);
	REQUIRE(allele_ids.size() == 1);
	REQUIRE(allele_ids[0] == 1);

	// insert allele 1 as an empty allele (should no longer contain previously inserted kmers)
	u.insert_empty_allele(1);
	REQUIRE(!u.kmer_on_path(0,1));
	REQUIRE(u.size() == 1);
	path_ids.clear();
	allele_ids.clear();
	u.get_path_ids(path_ids, allele_ids);
	REQUIRE(path_ids.size() == 1);
	REQUIRE(path_ids[0] == 1);
	REQUIRE(allele_ids.size() == 1);
	REQUIRE(allele_ids[0] == 1);
}

TEST_CASE("UniqueKmers kmers_on_alleles1", "[UniqueKmers kmers_on_alleles1]"){
	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1),CopyNumber(0.1,0.2,0.7)};
	UniqueKmers u(0,1000);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(cns[i], alleles[i]);
	}

	map<unsigned char, int> counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 2);
	REQUIRE(counts[1] == 2);
}


TEST_CASE("UniqueKmers kmers_on_alleles2", "[UniqueKmers kmers_on_alleles2]") {
	UniqueKmers u (0, 1000);
	u.insert_empty_allele(0);
	u.insert_empty_allele(2);
	u.insert_path(0,0);
	u.insert_path(1,0);
	map<unsigned char, int> counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 0);
	REQUIRE(counts[2] == 0);

	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1)};
	vector<vector<unsigned char>> alleles = { {2}, {0} };
	u.insert_kmer (cns[0], alleles[0]);
	counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 0);
	REQUIRE(counts[2] == 1);

	u.insert_kmer (cns[1], alleles[1]);
	counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 1);
	REQUIRE(counts[2] == 1);

	REQUIRE(!u.kmer_on_path(0,0));
	REQUIRE(!u.kmer_on_path(0,1));
	REQUIRE(u.kmer_on_path(1,0));
	REQUIRE(u.kmer_on_path(1,1));
}

TEST_CASE("UniqueKmers kmers_on_alleles3", "[UniqueKmers kmers_on_alleles3]") {
	UniqueKmers u(0, 1000);
	CopyNumber cn(0.9, 0.1, 0.2);
	vector<unsigned char> allele = {1};

	u.insert_kmer(cn, allele);
	u.insert_path(1,1);
	REQUIRE(u.size() == 1);
	REQUIRE(u.kmer_on_path(0, 1));
	map<unsigned char, int> counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 1);
	REQUIRE(counts[1] == 1);

	// insert allele 1 as an empty allele (should no longer contain previously inserted kmers)
	u.insert_empty_allele(1);
	REQUIRE(!u.kmer_on_path(0,1));
	counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 1);
	REQUIRE(counts[1] == 0);
}

TEST_CASE("UniqueKmers undefined_allele", "[UniqueKmers undefined_allele]"){
	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1),CopyNumber(0.1,0.2,0.7)};
	vector<vector<size_t>> paths = {{0,1,2}, {0,1}, {2}};
	UniqueKmers u(0,1000);

	// insert empty alleles
	u.insert_empty_allele(0, true);
	u.insert_empty_allele(1);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(cns[i], alleles[i]);
		for (auto& p : paths[i]){
			REQUIRE(u.kmer_on_path(i, p));
		}
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(u.get_variant_index() == 0);
		REQUIRE(u.get_variant_position() == 1000);
		REQUIRE(u.get_copynumber_of(i) == cns[i]);
	}

	REQUIRE(u.is_undefined_allele(0));
	REQUIRE(!u.is_undefined_allele(1));

	// get only defined alleles
	vector<unsigned char> all_alleles;
	vector<unsigned char> defined_alleles;
	u.get_allele_ids(all_alleles);
	REQUIRE(all_alleles.size() == 2);
	u.get_defined_allele_ids(defined_alleles);
	REQUIRE(defined_alleles.size() == 1);
	REQUIRE(defined_alleles.at(0) == (unsigned char) 1);

	// manually set allele to undefined
	u.set_undefined_allele(1);
	REQUIRE(u.is_undefined_allele(1));

	defined_alleles.clear();
	// now there are no defined alleles, so expect an empty vector
	u.get_defined_allele_ids(defined_alleles);
	REQUIRE(defined_alleles.empty());

	// make sure an execption is thrown in case an allele does not exist
	REQUIRE_THROWS(u.set_undefined_allele(2));
}
