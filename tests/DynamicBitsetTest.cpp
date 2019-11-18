#include "catch.hpp"
#include "../src/dynamicbitset.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("DynamicBitset empty", "[DynamicBitset empty]"){
	DynamicBitset b;
	REQUIRE(b.count() == 0);
	REQUIRE(b.is_set(10) == false);
	REQUIRE(b.is_set(65) == false);
	REQUIRE(b.convert_to_string() == "0000000000000000000000000000000000000000000000000000000000000000");
}

TEST_CASE("DynamicBitset set/unset", "[DynamicBitset set/unset]") {
	DynamicBitset b;
	b.set(0);
	REQUIRE(b.is_set(0));
	b.set(1);
	REQUIRE(b.is_set(1));
	REQUIRE(b.convert_to_string() == "1100000000000000000000000000000000000000000000000000000000000000"); 
	b.set(13);
	REQUIRE(b.is_set(13));
	REQUIRE(b.convert_to_string() == "1100000000000100000000000000000000000000000000000000000000000000");
	b.unset(13);
	REQUIRE(b.convert_to_string() == "1100000000000000000000000000000000000000000000000000000000000000");
	b.unset(14);
	REQUIRE(b.convert_to_string() == "1100000000000000000000000000000000000000000000000000000000000000");
	b.set(66);
	REQUIRE(b.is_set(66));
	b.unset(66);
	REQUIRE(!b.is_set(66));
	b.set(13);
	b.set(65);
	b.set(11);
}

TEST_CASE("DynamicBitset count", "DynamicBitset count") {
	DynamicBitset b;
	b.set(10);
	b.set(2);
	b.set(20);
	REQUIRE(b.count() == 3);
	b.unset(20);
	REQUIRE(b.count() == 2);
	b.set(1);
	REQUIRE(b.count() == 3);
	b.set(65);
	REQUIRE(b.count() == 4);
}

TEST_CASE("DynamicBitset operator&", "DynamicBitset operator&") {
	DynamicBitset a;
	DynamicBitset b;
	vector<unsigned int> bits_a = {0,2,4,6,8};
	vector<unsigned int> bits_b = {1,3,5,7,9};
	for (size_t i = 0; i < 5; ++i) {
		a.set(bits_a[i]);
		b.set(bits_b[i]);
	}	
	DynamicBitset result;
	result = a & b;
	REQUIRE(result.convert_to_string() == "0000000000000000000000000000000000000000000000000000000000000000");
	REQUIRE( result.count() == 0);
	result = b & b;
	REQUIRE(result.convert_to_string() == "0101010101000000000000000000000000000000000000000000000000000000");
}
