#include "catch.hpp"
#include "../src/uniquekmers.hpp"
#include "../src/copynumber.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("UniqueKmers testcase 1", "[UniqueKmers testcase 1]"){
	vector<string> kmers = {"ATGCA", "TTGAC", "GGTGA"};
	vector<CopyNumber> cns = {CopyNumber(0.05,0.9,0.05),CopyNumber(0.8,0.1,0.1),CopyNumber(0.1,0.2,0.7)};
	vector< vector<size_t> > paths = { {0,1,2}, {0,1}, {2} };
	UniqueKmers u(0, 1000);
	for (size_t i = 0; i < 3; ++i){
		u.insert_kmer(kmers[i], cns[i], paths[i]);
		for (auto& p : paths[i]){
			REQUIRE(u.kmer_on_path(kmers[i], p));
		}
	}
	for (size_t i = 0; i < 3; ++i){
		REQUIRE(u.get_variant_index() == 0);
		REQUIRE(u.get_variant_position() == 1000);
		REQUIRE(u.contains_kmer(kmers[i]));
		REQUIRE(u.get_copynumber_of(kmers[i]) == cns[i]);
	}
}

TEST_CASE("UniqueKmers get_copynumber_of", "[UniqueKmers get_copynumber_of]"){
	UniqueKmers u(1, 2000);
	vector<string> kmers = {"ATGCA", "TTGAC", "GGTGA"};
	for (size_t i = 0; i < 3; ++i){
		CHECK_THROWS(u.get_copynumber_of(kmers[i]));
	}
}

TEST_CASE("UniqueKmers insert_empty_path", "[UniqueKmers insert_empty_path]") {
	UniqueKmers u(0, 1000);
	vector<string> kmers = {"ATAG", "CTCG", "GGGG"};
	vector<CopyNumber> cns = { CopyNumber(0.001, 0.5, 0.001), CopyNumber(1.0, 0.0, 0.0), CopyNumber(0.001,0.001, 0.6) };
	vector<vector<size_t>> paths = { {0}, {0}, {2} };

	// insert the kmers
	for (size_t i = 0; i < 3; ++i) {
		u.insert_kmer(kmers[i], cns[i], paths[i]);
		for (auto& p : paths[i]) {
			REQUIRE(u.kmer_on_path(kmers[i], p));
		}
	}

	vector<size_t> path_ids;
	u.get_path_ids(path_ids);
	REQUIRE(u.size() == 3);

	REQUIRE(path_ids.size() == 2);
	REQUIRE(path_ids[0] == 0);
	REQUIRE(path_ids[1] == 2);

	// insert an empty path
	u.insert_empty_path(1);
	path_ids.clear();
	u.get_path_ids(path_ids);
	REQUIRE(u.size() == 3);

	REQUIRE(path_ids.size() == 3);
	REQUIRE(path_ids[0] == 0);
	REQUIRE(path_ids[1] == 1);
	REQUIRE(path_ids[2] == 2);

	// make sure path is indeed empty
	for (size_t i = 0; i < 3; ++i) {
		REQUIRE(! u.kmer_on_path(kmers[i], 1));
	}

	// add kmer to the empty path
	vector<size_t> path = {1};
	u.insert_kmer("CCCC", CopyNumber(0.01, 0.8, 0.01), path);
	REQUIRE(u.size() == 4);
	REQUIRE(u.kmer_on_path("CCCC", 1));
}

TEST_CASE("UniqueKmers insert_empty_path2", "[UniqueKmers insert_empty_path2]") {
	UniqueKmers u(0, 1000);
	string kmer = "AAAAA";
	CopyNumber cn(0.9, 0.1, 0.2);
	vector<size_t> path = {1};

	u.insert_kmer(kmer, cn, path);
	REQUIRE(u.size() == 1);
	REQUIRE(u.kmer_on_path(kmer, 1));

	vector<size_t> path_ids;
	u.get_path_ids(path_ids);
	REQUIRE(path_ids.size() == 1);
	REQUIRE(path_ids[0] == 1);

	// insert path 1 as an empty path (should no longer contain previously inserted path)
	u.insert_empty_path(1);
	REQUIRE(!u.kmer_on_path(kmer,1));
	REQUIRE(u.size() == 1);
	path_ids.clear();
	u.get_path_ids(path_ids);
	REQUIRE(path_ids.size() == 1);
	REQUIRE(path_ids[0] == 1);

}
