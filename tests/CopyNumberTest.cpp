#include "catch.hpp"
#include "utils.hpp"
#include "../src/copynumber.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("CopyNumber get_probability_of", "[CopyNumber get_probability_of]"){
	CopyNumber c(0.9,0.1,0.0);
	REQUIRE(doubles_equal(c.get_probability_of(0), 0.9));
	CopyNumber c2(0.4,0.5,0.1);
	REQUIRE(doubles_equal(c2.get_probability_of(1), 0.5));
	CopyNumber c3(0.0,0.0,1.0);
	REQUIRE(doubles_equal(c3.get_probability_of(2), 1.0));
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

TEST_CASE("CopyNumber scaling",  "[CopyNumber scaling]") {
	// regularization const. of 0.0 leads to normalizing the probabilities to sum up to 1
	CopyNumber c1(0.1, 0.1, 0.1, 0.0);
	REQUIRE(doubles_equal(c1.get_probability_of(0), 1.0/3.0));
	REQUIRE(doubles_equal(c1.get_probability_of(1), 1.0/3.0));
	REQUIRE(doubles_equal(c1.get_probability_of(2), 1.0/3.0));

	CopyNumber c2(0.001, 0.6, 0.0004, 0.0);
	REQUIRE(doubles_equal(c2.get_probability_of(0), 0.001 / 0.6014));
	REQUIRE(doubles_equal(c2.get_probability_of(1), 0.6 / 0.6014));
	REQUIRE(doubles_equal(c2.get_probability_of(2), 0.0004 / 0.6014));
}

TEST_CASE("CopyNumber regularization", "[CopyNumber regularization]") {
	CopyNumber c1(0.1, 0.1, 0.1, 10);
	REQUIRE(doubles_equal(c1.get_probability_of(0), 1.0/3.0));
	REQUIRE(doubles_equal(c1.get_probability_of(1), 1.0/3.0));
	REQUIRE(doubles_equal(c1.get_probability_of(2), 1.0/3.0));

	CopyNumber c2(0.2, 0.9, 1.1, 100.0);
	REQUIRE(doubles_equal(c2.get_probability_of(0), 0.33156849768));
	REQUIRE(doubles_equal(c2.get_probability_of(1), 0.33388484447));
	REQUIRE(doubles_equal(c2.get_probability_of(2), 0.33454665784));

	CopyNumber c3(0.001, 0.0, 0.999, 20);
	REQUIRE(doubles_equal(c3.get_probability_of(0), 0.3278852459));
	REQUIRE(doubles_equal(c3.get_probability_of(1), 0.32786885245));
	REQUIRE(doubles_equal(c3.get_probability_of(2), 0.34424590163));
	
}


