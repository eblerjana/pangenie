#include <math.h>
#include "catch.hpp"
#include "utils.hpp"
#include "../src/collapsedtransitionprobabilitycomputer.hpp"
#include "../src/uniquekmers.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("CollapsedTransitionProbabilityComputer compute_transition_prob1", "[TransitionProbabilityComputer compute_transition_probability]") {
	CollapsedTransitionProbabilityComputer t (1000000, 2000000, 1.26, 5, false, 0.25);
	double recomb_prob =  0.04455105238;
	double no_recomb_prob = recomb_prob + 0.77724473806;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	// first variant: h1 - 0, h2+h3 - 1, second variant: h1,h2 - 0, h3 - 1
	UniqueKmers u1 (0, 1000000);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_path(0,0);
	u1.insert_path(1,1);
	u1.insert_path(2,1);

	UniqueKmers u2 (1, 2000000);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);
	u2.insert_path(0,0);
	u2.insert_path(1,0);
	u2.insert_path(2,1);

	// paths for genotypes 0|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*2 + two_recomb*1));
	// paths for genotypes 0|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*4 + two_recomb*3));
	// paths for genotypes 1|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*4 + two_recomb*3));
	// paths for genotypes 1|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*6 + two_recomb*9));

	// paths for genotypes 0|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*1));
	// paths for genotypes 0|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*1 + one_recomb*2 + two_recomb*1));
	// paths for genotypes 1|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*3));
	// paths for genotypes 1|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*1 + one_recomb*4 + two_recomb*3));

	// paths for genotypes 0|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*1));
	// paths for genotypes 0|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*3));
	// paths for genotypes 1|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*2 + two_recomb*1));
	// paths for genotypes 1|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*4 + two_recomb*3));

	// paths for genotypes 0|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*1));
	// paths for genotypes 0|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*1));
	// paths for genotypes 1|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*1));
	// paths for genotypes 1|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*1 + one_recomb*2 + two_recomb*1));
}

TEST_CASE("CollapsedTransitionProbabilityComputer compute_transition_prob2", "[CollapsedTransitionProbabilityComputer compute_transition_probability2]") {
	CollapsedTransitionProbabilityComputer t (1000000, 2000000, 1.26, 10, false, 0.25);
	double recomb_prob = 0.01183851532;
	double no_recomb_prob = recomb_prob + 0.88161484678;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	// first variant: h1 - 0, h2 - 1, second variant: h2 - 0, h1 - 1
	UniqueKmers u1 (0, 1000000);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_path(0,0);
	u1.insert_path(1,1);

	UniqueKmers u2 (1, 2000000);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);
	u2.insert_path(0,1);
	u2.insert_path(1,0);


	// paths for genotypes 0|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*1));
	// paths for genotypes 0|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 1|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 1|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*0 + two_recomb*0));

	// paths for genotypes 0|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 0|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*1));
	// paths for genotypes 1|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*1 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*0));

	// paths for genotypes 0|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 0|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*1 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*1));
	// paths for genotypes 1|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*1 + two_recomb*0));

	// paths for genotypes 0|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*1 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 0|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 1|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*1 + two_recomb*0));
	// paths for genotypes 1|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*1));
}


TEST_CASE("CollapsedTransitionProbabilityComputer compute_transition_prob3", "[CollapsedTransitionProbabilityComputer compute_transition_probability3]") {
	CollapsedTransitionProbabilityComputer t (1000000, 2000000, 1.26, 10, false, 0.25);
	double recomb_prob = 0.01183851532;
	double no_recomb_prob = recomb_prob + 0.88161484678;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	// first variant: h1 - 0, second variant: h2 - 1
	UniqueKmers u1 (0, 1000000);
	u1.insert_empty_allele(0);
	u1.insert_empty_allele(1);
	u1.insert_path(0,0);

	UniqueKmers u2 (1, 2000000);
	u2.insert_empty_allele(0);
	u2.insert_empty_allele(1);
	u2.insert_path(1,1);

	// paths for genotypes 0|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 0|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|0 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|1 -> 0|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));

	// paths for genotypes 0|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 0|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|0 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|1 -> 0|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(0), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));

	// paths for genotypes 0|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 0|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|0 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|1 -> 1|0
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(0)), no_recomb*0 + one_recomb*0 + two_recomb*0));

	// paths for genotypes 0|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*1));
	// paths for genotypes 0|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(0), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|0 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(0), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));
	// paths for genotypes 1|1 -> 1|1
	REQUIRE(doubles_equal(t.compute_transition_prob(u1.get_paths_of_allele(1), u1.get_paths_of_allele(1), u2.get_paths_of_allele(1), u2.get_paths_of_allele(1)), no_recomb*0 + one_recomb*0 + two_recomb*0));

}

TEST_CASE("CollapsedTransitionProbabilityComputer compute_transition_start", "[CollapsedTransitioProbabilityComputer compute_transition_start]") {
	CollapsedTransitionProbabilityComputer t (1000000, 2000000, 1.26, 10, false, 0.25);
	double recomb_prob = 0.01183851532;
	double no_recomb_prob = recomb_prob + 0.88161484678;

	double no_recomb = no_recomb_prob * no_recomb_prob;
	double one_recomb = recomb_prob * no_recomb_prob;
	double two_recomb = recomb_prob*recomb_prob;

	UniqueKmers u (0,1000000);
	u.insert_path(0,0);
	u.insert_path(1,0);
	u.insert_path(2,1);

	// transition probabilities are not normalized (as normalization factor is a constant and cancels out in HMM)
	REQUIRE(doubles_equal(t.compute_transition_start(u.get_paths_of_allele(0), u.get_paths_of_allele(0)), 4.0));
	REQUIRE(doubles_equal(t.compute_transition_start(u.get_paths_of_allele(0), u.get_paths_of_allele(1)), 2.0));
	REQUIRE(doubles_equal(t.compute_transition_start(u.get_paths_of_allele(1), u.get_paths_of_allele(0)), 2.0));
	REQUIRE(doubles_equal(t.compute_transition_start(u.get_paths_of_allele(1), u.get_paths_of_allele(1)), 1.0));
}
