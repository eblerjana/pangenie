#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/probabilitycomputer.hpp"

#include <iostream>

using namespace std;

TEST_CASE ("ProbabilityComputer get_probability", "[ProbabilityComputer get_probability]") {
	ProbabilityComputer p(0.1, 10, 50);
	REQUIRE(doubles_equal(p.get_probability(2,70), 0.00136386));
	REQUIRE(doubles_equal(p.get_probability(0,2), 0.081));
	REQUIRE(doubles_equal(p.get_probability(1,10), 0.12511));
}

TEST_CASE ("ProbabilityComputer set_parameters", "[ProbabilityComputer set_parameters]") {
	ProbabilityComputer p;
	// parameters not yet set
	CHECK_THROWS (p.get_probability(1, 20));
	p.set_parameters (1, 25, 50);
	REQUIRE (doubles_equal(p.get_probability(1,20), 0.0519175));
}
