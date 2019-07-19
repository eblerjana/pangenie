#include "catch.hpp"
#include "../src/columnindexer.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("ColumnIndexer insert", "[ColumnIndexer insert]"){
	ColumnIndexer c(0, 3);
	c.insert(make_pair(0,0), make_pair(0,0));
	REQUIRE(c.get_size() == 1);
	c.insert(make_pair(0,1), make_pair(0,2));
	REQUIRE(c.get_size() == 2);
	c.insert(make_pair(1,0), make_pair(2,0));
	REQUIRE(c.get_size() == 3);
	REQUIRE(c.get_paths(0) == pair<size_t,size_t>(0,0));
	REQUIRE(c.get_paths(1) == pair<size_t,size_t>(0,1));
	REQUIRE(c.get_paths(2) == pair<size_t,size_t>(1,0));
	REQUIRE(c.get_alleles(0) == pair<unsigned char, unsigned char>(0,0));
	REQUIRE(c.get_alleles(1) == pair<unsigned char, unsigned char>(0,2));
	REQUIRE(c.get_alleles(2) == pair<unsigned char, unsigned char>(2,0));
}

TEST_CASE("ColumnIndexer test_invalid", "[ColumnIndexer test_invalid]") {
	ColumnIndexer c(0, 3);
	REQUIRE(c.get_size() == 0);
	REQUIRE_THROWS(c.get_paths(0));
	REQUIRE_THROWS(c.get_paths(1));
	REQUIRE_THROWS(c.get_alleles(0));
	REQUIRE_THROWS(c.get_alleles(1));

	c.insert(make_pair(1,2), make_pair(0,0));
	REQUIRE(c.get_size() == 1);
	REQUIRE(c.get_paths(0) == pair<size_t,size_t>(1,2));
	REQUIRE(c.get_alleles(0) == pair<unsigned char, unsigned char>(0,0));
}
