#include "catch.hpp"
#include "../src/columnindexer.hpp"
#include "../src/uniquekmers.hpp"
#include <vector>
#include <string>
#include <memory>

using namespace std;

TEST_CASE("ColumnIndexer testcase 1", "[ColumnIndexer testcase 1]"){
	vector<unsigned char> path_to_allele = {0, 1, 0, 0, 0};
	shared_ptr<UniqueKmers> u1 = shared_ptr<UniqueKmers>(new UniqueKmers (2000, path_to_allele));
	vector<unsigned char> a1 = {0};
	vector<unsigned char> a2 = {1};
	u1->insert_kmer(10, a1);
	u1->insert_kmer(10, a2);
	u1->set_coverage(5);

	// position only covered by reference alleles and should be skipped
	path_to_allele = {0, 0, 1, 1, 1};
	shared_ptr<UniqueKmers> u2 =  shared_ptr<UniqueKmers> ( new UniqueKmers (2500, path_to_allele));
	u2->insert_kmer(10, a1);
	u2->insert_kmer(20, a2);

	path_to_allele = {0, 0, 1, 1, 1};
	shared_ptr<UniqueKmers> u3 = shared_ptr<UniqueKmers> (new UniqueKmers (3000, path_to_allele));
	u3->insert_kmer(20, a1);
	u3->insert_kmer(5, a2);
	u3->set_coverage(5);

	vector<shared_ptr<UniqueKmers>> unique_kmers = {u1, u2, u3};
	vector<unsigned short> only_paths = {2,3};

	ColumnIndexer column_indexer(&unique_kmers, &only_paths);
	// first variant only covered by reference paths
	REQUIRE(column_indexer.size() == 2);
	REQUIRE(column_indexer.get_variant_id(0) == 1);
	REQUIRE(column_indexer.get_variant_id(1) == 2);
	REQUIRE(column_indexer.nr_paths() == 2);
	REQUIRE(column_indexer.get_path(0) == 2);
	REQUIRE(column_indexer.get_path(1) == 3);
	REQUIRE(column_indexer.get_allele(0,0) == 1);
	REQUIRE(column_indexer.get_allele(1,0) == 1);
	REQUIRE(column_indexer.get_allele(0,1) == 1);
	REQUIRE(column_indexer.get_allele(1,1) == 1);

	column_indexer = ColumnIndexer(&unique_kmers, nullptr);
	REQUIRE(column_indexer.size() == 3);
	REQUIRE(column_indexer.get_variant_id(0) == 0);
	REQUIRE(column_indexer.get_variant_id(1) == 1);
	REQUIRE(column_indexer.get_variant_id(2) == 2);
	REQUIRE(column_indexer.nr_paths() == 5);
	REQUIRE(column_indexer.get_path(0) == 0);
	REQUIRE(column_indexer.get_path(1) == 1);
	REQUIRE(column_indexer.get_path(2) == 2);
	REQUIRE(column_indexer.get_path(3) == 3);
	REQUIRE(column_indexer.get_path(4) == 4);
	REQUIRE(column_indexer.get_allele(0,0) == 0);
	REQUIRE(column_indexer.get_allele(1,0) == 1);
	REQUIRE(column_indexer.get_allele(2,1) == 1);
	REQUIRE(column_indexer.get_allele(4,2) == 1);
	REQUIRE(column_indexer.get_allele(1,2) == 0);

	CHECK_THROWS(column_indexer.get_variant_id(3));
	CHECK_THROWS(column_indexer.get_path(5));
	CHECK_THROWS(column_indexer.get_allele(3,3));
	CHECK_THROWS(column_indexer.get_allele(5,1));
}