#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/samplingtransitions.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("SamplingTransitions compute_transition_cost1", "[SamplingTransitions compute_transition_cost1]") {

	SamplingTransitions s (1000000, 2000000, 1.26, 5, 0.25);
	double recomb_prob =  0.04455105238;
	double no_recomb_prob = recomb_prob + 0.77724473806;

	unsigned int expected_cost = -10.0 * log10(recomb_prob);

	REQUIRE(s.compute_transition_cost(false) == 0);
	REQUIRE(s.compute_transition_cost(true) == expected_cost);
}
