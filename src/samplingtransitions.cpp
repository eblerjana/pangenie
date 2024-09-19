#include <cmath>
#include <cassert>
#include "samplingtransitions.hpp"

SamplingTransitions::SamplingTransitions(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, long double effective_N) {
	assert(from_variant <= to_variant);
	// using same formula as in WhatsHap
//	long double distance = (to_variant - from_variant) * 0.000001 * ((long double) recomb_rate) * 4.0L * effective_N;
	long double distance = (to_variant - from_variant) * 0.000004L * ((long double) recomb_rate) * effective_N;
	// use Li-Stephans pair HMM transitions TODO: correct?
	long double recomb_prob = (1.0L - exp(-distance / (long double) nr_paths) )* (1.0L / (long double) nr_paths);
	this->cost = -10.0 * log10(recomb_prob);
}

unsigned int SamplingTransitions::compute_transition_cost(bool recombination) {
	if (recombination) {
		return this->cost;
	} else {
		return 0;
	}
}
