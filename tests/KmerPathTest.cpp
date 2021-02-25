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
	REQUIRE(p1.convert_to_string() == "00000000000000000000010000000000");
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

TEST_CASE("KmerPath set_position", "[KmerPath set_position]") {
	KmerPath p1;
	p1.set_position(0);
	REQUIRE(p1.convert_to_string() == "10000000000000000000000000000000");
	p1.set_position(32);
	REQUIRE(p1.convert_to_string() == "1000000000000000000000000000000010000000000000000000000000000000");
}
