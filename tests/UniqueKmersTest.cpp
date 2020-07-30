#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("UniqueKmers testcase 1", "[UniqueKmers testcase 1]"){
	vector<unsigned short> counts = {5,6,7};
	vector<vector<size_t>> paths = {{0,1,2}, {0,1}, {2}};
	UniqueKmers u(1000);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(counts[i], alleles[i]);
		for (auto& p : paths[i]){
			REQUIRE(u.kmer_on_path(i, p));
		}
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(u.get_variant_position() == 1000);
		REQUIRE(u.get_readcount_of(i) == counts[i]);
	}
}

TEST_CASE("UniqueKmers get_copynumber_of", "[UniqueKmers get_copynumber_of]"){
	UniqueKmers u(2000);
	for (size_t i = 0; i < 3; ++i){
		CHECK_THROWS(u.get_readcount_of(i));
	}
	vector<unsigned char> alleles = {0};
	u.insert_kmer(5, alleles);
	REQUIRE(u.get_readcount_of(0) == 5);
}

TEST_CASE("UniqueKmers insert_empty_path", "[UniqueKmers insert_empty_path]") {
	UniqueKmers u(1000);
	vector<unsigned short> counts = {5, 0, 10};
	vector<vector<size_t>> paths = { {0}, {0}, {2} };
	vector<vector<unsigned char>> alleles = {{0}, {0}, {1}};

	u.insert_path(0,0);
	u.insert_path(2,1);

	// insert the kmers
	for (size_t i = 0; i < 3; ++i) {
		u.insert_kmer(counts[i], alleles[i]);
		for (auto& p : paths[i]) {
			REQUIRE(u.kmer_on_path(i, p));
		}
	}

	vector<unsigned short> path_ids;
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
	for (unsigned short i = 0; i < 3; ++i) {
		REQUIRE(! u.kmer_on_path(i, 1));
	}

	// add kmer to the empty path
	vector<unsigned char> allele = {2};
	u.insert_kmer(5, allele);
	REQUIRE(u.size() == 4);
	REQUIRE(u.kmer_on_path(3, 1));
}

TEST_CASE("UniqueKmers insert_empty_path2", "[UniqueKmers insert_empty_path2]") {
	UniqueKmers u(1000);
	vector<unsigned char> allele = {1};

	u.insert_kmer(10, allele);
	u.insert_path(1,1);
	REQUIRE(u.size() == 1);
	REQUIRE(u.kmer_on_path(0, 1));

	vector<unsigned short> path_ids;
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
	vector<unsigned short> read_counts = {5, 1, 9};
	UniqueKmers u(1000);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(read_counts[i], alleles[i]);
	}

	map<unsigned char, int> counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 2);
	REQUIRE(counts[1] == 2);
}


TEST_CASE("UniqueKmers kmers_on_alleles2", "[UniqueKmers kmers_on_alleles2]") {
	UniqueKmers u (1000);
	u.insert_empty_allele(0);
	u.insert_empty_allele(2);
	u.insert_path(0,0);
	u.insert_path(1,0);
	map<unsigned char, int> counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 0);
	REQUIRE(counts[2] == 0);

	vector<unsigned short> read_counts = {5, 1};
	vector<vector<unsigned char>> alleles = { {2}, {0} };
	u.insert_kmer (read_counts[0], alleles[0]);
	counts = u.kmers_on_alleles();
	REQUIRE(counts.size() == 2);
	REQUIRE(counts[0] == 0);
	REQUIRE(counts[2] == 1);

	u.insert_kmer (read_counts[1], alleles[1]);
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
	UniqueKmers u(1000);
	vector<unsigned char> allele = {1};

	u.insert_kmer(10, allele);
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

TEST_CASE("UniqueKmers get_path_ids", "[UniqueKmers get_path_ids]") {
	UniqueKmers u (1000);
	u.insert_empty_allele(0);
	u.insert_empty_allele(2);
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,2);
	u.insert_path(3,1);

	vector<unsigned short> path_ids;
	vector<unsigned char> allele_ids;
	u.get_path_ids(path_ids, allele_ids);

	vector<unsigned short> expected_path_ids = {0,1,2,3};
	vector<unsigned char> expected_allele_ids = {0,0,2,1};
	REQUIRE(path_ids == expected_path_ids);
	REQUIRE(allele_ids == expected_allele_ids);

	// select only specific path_ids
	path_ids.clear();
	allele_ids.clear();
	vector<unsigned short> specific_ids = {0,2,10};
	expected_path_ids = {0,2};
	expected_allele_ids = {0,2};
	u.get_path_ids(path_ids, allele_ids, &specific_ids);
	REQUIRE(path_ids == expected_path_ids);
	REQUIRE(allele_ids == expected_allele_ids);

	// no overlap between requested ids and the ones in UniqueKmers
	path_ids.clear();
	allele_ids.clear();
	specific_ids = {20,30,40};
	expected_path_ids = {};
	expected_allele_ids = {};
	u.get_path_ids(path_ids, allele_ids, &specific_ids);
	REQUIRE(path_ids == expected_path_ids);
	REQUIRE(allele_ids == expected_allele_ids);

	// all alleles should have been returned
	path_ids.clear();
	allele_ids.clear();
	specific_ids = {0,1,2,3};
	expected_path_ids = {0,1,2,3};
	expected_allele_ids = {0,0,2,1};
	u.get_path_ids(path_ids, allele_ids, &specific_ids);
	REQUIRE(path_ids == expected_path_ids);
	REQUIRE(allele_ids == expected_allele_ids);	
}

TEST_CASE("UniqueKmers undefined_allele", "[UniqueKmers undefined_allele]"){
	vector<unsigned short> read_counts = {10, 1, 20};
	vector<vector<size_t>> paths = {{0,1,2}, {0,1}, {2}};
	UniqueKmers u(1000);

	// insert empty alleles
	u.insert_empty_allele(0, true);
	u.insert_empty_allele(1);

	// insert paths
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	vector< vector<unsigned char> > alleles = { {0,1}, {0}, {1} };
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(read_counts[i], alleles[i]);
		for (auto& p : paths[i]){
			REQUIRE(u.kmer_on_path(i, p));
		}
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(u.get_variant_position() == 1000);
		REQUIRE(u.get_readcount_of(i) == read_counts[i]);
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
