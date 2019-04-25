#include "catch.hpp"
#include "../src/copynumber.h"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("CopyNumber simple example", "[CopyNumber simple]"){
	CopyNumber c(0.9,0.1,0.0);
	REQUIRE(c.get_copynumber() == 0);
	CopyNumber c2(0.4,0.5,0.1);
	REQUIRE(c2.get_copynumber() == 1);
	CopyNumber c3(0.0,0.0,1.0);
	REQUIRE(c3.get_copynumber() == 2);
}

TEST_CASE("CopyNumber operators", "[CopyNumber operators]"){
	CopyNumber c1(0.1,0.2,0.7);
	CopyNumber c2(0.1,0.2,0.7);
	CopyNumber c3(0.0,1.0,0.0);

	REQUIRE(c1 == c2);
	REQUIRE(c2 == c1);
	REQUIRE_FALSE(c1 == c3);
	REQUIRE_FALSE(c2 == c3);
	REQUIRE(c1 != c3);
	REQUIRE(c2 != c3);
}
