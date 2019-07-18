#include <cassert>
#include <math.h>
#include "transitionprobabilitycomputer.hpp"
#include <iostream>

using namespace std;

TransitionProbabilityComputer::TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate)
{
	assert(from_variant <= to_variant );
	// using same formula as in WhatsHap
	long double distance = (to_variant - from_variant) * 0.000001 * (long double) recomb_rate;
	// using formula from: https://en.wikipedia.org/wiki/Centimorgan
	long double recomb_prob = (1.0L - exp(-(2.0L*distance)/100)) / 2.0L;
	long double no_recomb_prob = 1.0L - recomb_prob;
	this->probabilities = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};
}

long double TransitionProbabilityComputer::compute_transition_prob(size_t path_id1, size_t path_id2, size_t path_id3, size_t path_id4){
	// determine number of recombination events
	unsigned int nr_events = 0;
	if (path_id1 != path_id3) nr_events += 1;
	if (path_id2 != path_id4) nr_events += 1;
	return this->probabilities[nr_events];
}
