#include "catch.hpp"
#include "../src/columnindexer.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("ColumnIndexer insert", "[ColumnIndexer insert]"){
	ColumnIndexer c(1);
	c.insert_path(0,0);
	c.insert_path(3,2);
	REQUIRE(c.nr_paths() == 2);
	REQUIRE(c.get_path(0) == 0);
	REQUIRE(c.get_path(1) == 3);
	REQUIRE(c.get_allele(0) == 0);
	REQUIRE(c.get_allele(1) == 2);
	REQUIRE(c.get_path_ids_at(0) == pair<unsigned short,unsigned short>(0,0));
	REQUIRE(c.get_path_ids_at(1) == pair<unsigned short,unsigned short>(0,1));
	REQUIRE(c.get_path_ids_at(2) == pair<unsigned short,unsigned short>(1,0));
	REQUIRE(c.get_path_ids_at(3) == pair<unsigned short,unsigned short>(1,1));
	REQUIRE(c.get_variant_id() == 1);
}

TEST_CASE("ColumnIndexer test_invalid", "[ColumnIndexer test_invalid]") {
	ColumnIndexer c(0);
	REQUIRE(c.nr_paths() == 0);
	REQUIRE_THROWS(c.get_path(0));
	REQUIRE_THROWS(c.get_path(1));
	REQUIRE_THROWS(c.get_allele(0));
	REQUIRE_THROWS(c.get_allele(1));
	REQUIRE_THROWS(c.get_path_ids_at(0));

	c.insert_path(1,0);
	REQUIRE(c.nr_paths() == 1);
	REQUIRE(c.get_path(0) == 1);
	REQUIRE(c.get_allele(0) == 0);
	REQUIRE(c.get_path_ids_at(0) == pair<unsigned short,unsigned short>(0,0));
	REQUIRE(c.get_variant_id() == 0);
}
