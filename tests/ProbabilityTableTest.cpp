#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/probabilitytable.hpp"

#include <iostream>

using namespace std;

TEST_CASE ("ProbabilityTable1", "[ProbabilityTable1]") {
	ProbabilityTable p(5,6,1,0.0);

	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(0), 0.99));
	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(1), 0.08208499862));
	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(2), 0.00673794699));

	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(0), 0.0099));
	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(1), 0.20521249655));
	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(2), 0.03368973499));

	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(0), 0.99));
	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(1), 0.04978706836));
	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(2), 0.00247875217));

	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(0), 0.0099));
	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(1), 0.149361205103));
	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(2), 0.014872513059));
}

TEST_CASE ("ProbabilityTable2", "[ProbabilityTable2]") {
	ProbabilityTable p(4,7,2,0.0);

	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(0), 0.99));
	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(1), 0.08208499862));
	REQUIRE(doubles_equal(p.get_probability(5,0).get_probability_of(2), 0.00673794699));

	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(0), 0.0099));
	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(1), 0.20521249655));
	REQUIRE(doubles_equal(p.get_probability(5,1).get_probability_of(2), 0.03368973499));

	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(0), 0.99));
	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(1), 0.04978706836));
	REQUIRE(doubles_equal(p.get_probability(6,0).get_probability_of(2), 0.00247875217));

	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(0), 0.0099));
	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(1), 0.149361205103));
	REQUIRE(doubles_equal(p.get_probability(6,1).get_probability_of(2), 0.014872513059));
}
