#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/transitionprobabilitycomputer.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("TransitionProbabilityComputer compute_transition_prob1", "[TransitionProbabilityComputer compute_transition_probability]") {
	TransitionProbabilityComputer t (1000000, 2000000, 1.26, 5, false, 0.25);
	double recomb_prob =  0.04455105238;
	double no_recomb_prob = recomb_prob + 0.77724473806;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,1), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,0,1,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,2,2,1), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,0), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,1), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(3,3,0,0), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,3,1,1), one_recomb));
}

TEST_CASE("TransitionProbabilityComputer compute_transition_prob2", "[TransitionProbabilityComputer compute_transition_probability]") {
	TransitionProbabilityComputer t (1000000, 2000000, 1.26, 10, false, 0.25);
	double recomb_prob = 0.01183851532;
	double no_recomb_prob = recomb_prob + 0.88161484678;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,0,1), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,0,1,0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,2,2,1), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,0), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(0,0,1,1), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(3,3,0,0), two_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1,3,1,1), one_recomb));
}

TEST_CASE("TransitionProbabilityComputer compute_transition_prob1_switches", "[TransitionProbabilityComputer compute_transition_probability_switches]") {
	TransitionProbabilityComputer t (1000000, 2000000, 1.26, 5, false, 0.25);
	double recomb_prob =  0.04455105238;
	double no_recomb_prob = recomb_prob + 0.77724473806;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	REQUIRE(doubles_equal(t.compute_transition_prob(0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(2), two_recomb));
}

TEST_CASE("TransitionProbabilityComputer compute_transition_prob2_switches", "[TransitionProbabilityComputer compute_transition_probability_switches]") {
	TransitionProbabilityComputer t (1000000, 2000000, 1.26, 10, false, 0.25);
	double recomb_prob = 0.01183851532;
	double no_recomb_prob = recomb_prob + 0.88161484678;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	REQUIRE(doubles_equal(t.compute_transition_prob(0), no_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(1), one_recomb));
	REQUIRE(doubles_equal(t.compute_transition_prob(2), two_recomb));
}