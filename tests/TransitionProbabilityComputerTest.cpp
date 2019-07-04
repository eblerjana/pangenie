#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/transitionprobabilitycomputer.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("TransitionProbabilityComputer compute_transition_probability", "[TransitionProbabilityComputer compute_transition_probability]") {
	TransitionProbabilityComputer t (1000000, 2000000, 1.26);
	float no_recomb = (1-0.01244256522) * (1-0.01244256522);
	float one_recomb = 0.01244256522 * (1 - 0.01244256522);
	float two_recomb = 0.01244256522 * 0.01244256522;

	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,1), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,0,1,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,2,2,1), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,0), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,1), two_recomb));
}
