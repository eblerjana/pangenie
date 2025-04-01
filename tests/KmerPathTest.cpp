#include "catch.hpp"
#include "../src/kmerpath.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("KmerPath getter/setter", "[KmerPath getter/setter]"){
	vector<string> kmers = {"ATGC", "AAAA", "GCTA"};
	KmerPath p1;
	KmerPath p2;
	for (size_t i = 0; i < kmers.size(); ++i){
		p1.set_position(i);
		REQUIRE(p1.get_position(i) == 1);
	}
	p2.set_position(0);
	REQUIRE(p2.get_position(0) == 1);
	REQUIRE(p1.convert_to_string() == "11100000000000000000000000000000");
	REQUIRE(p2.convert_to_string() == "10000000000000000000000000000000");
}

TEST_CASE("KmerPath getter/setter2", "[KmerPath getter/setter2]"){
	KmerPath p1;
	KmerPath p2;
	p1.set_position(21);
	REQUIRE(p1.get_position(21) == 1);
	p2.set_position(0);
	REQUIRE(p2.get_position(0) == 1);
	REQUIRE(p1.convert_to_string() == "00000000000000000000010000000000000000000000000000000");
	REQUIRE(p2.convert_to_string() == "10000000000000000000000000000000");
}

TEST_CASE("KmerPath set multiple", "[KmerPath set multiple]"){
	KmerPath p1;
	p1.set_position(0);
	REQUIRE(p1.convert_to_string() == "10000000000000000000000000000000");
	REQUIRE(p1.get_position(0) == 1);
	p1.set_position(0);
	REQUIRE(p1.convert_to_string() == "10000000000000000000000000000000");
	REQUIRE(p1.get_position(0) == 1);
}

TEST_CASE("KmerPath set_position_invalid", "[KmerPath set_position_invalid]") {
	KmerPath p1;
	p1.set_position(0);
	REQUIRE(p1.convert_to_string() == "10000000000000000000000000000000");
	REQUIRE_THROWS(p1.set_position(32));
}

TEST_CASE("KmerPath set_position_offset", "[KmerPath set_position_offset]") {
	KmerPath p;
	p.set_position(32);
	p.set_position(34);
	REQUIRE(p.get_position(32) == 1);
	REQUIRE(p.get_position(34) == 1);
	REQUIRE(p.convert_to_string() == "0000000000000000000000000000000010100000000000000000000000000000");
	CHECK_THROWS(p.set_position(31));
	REQUIRE(p.get_position(31) == 0);
	CHECK_THROWS(p.set_position(64));
	REQUIRE(p.get_position(64) == 0);
}

TEST_CASE("KmerPath set_position_offset2", "[KmerPath set_position_offset2]") {
	KmerPath p;
	p.set_position(31);
	p.set_position(32);
	REQUIRE(p.get_position(31) == 1);
	REQUIRE(p.get_position(32) == 1);
	REQUIRE(p.convert_to_string() == "000000000000000000000000000000011000000000000000000000000000000");
	CHECK_THROWS(p.set_position(64));
	REQUIRE(p.get_position(64) == 0);
}