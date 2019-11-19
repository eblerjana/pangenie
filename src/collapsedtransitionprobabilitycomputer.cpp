#include <cassert>
#include <math.h>
#include "collapsedtransitionprobabilitycomputer.hpp"
#include <iostream>

using namespace std;

CollapsedTransitionProbabilityComputer::CollapsedTransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, size_t nr_paths, bool uniform, long double effective_N)
{
	assert(from_variant <= to_variant );

	if (uniform) {
		this->probabilities = {1.0L, 1.0L, 1.0L, 1.0L};
	} else {
		// using same formula as in WhatsHap
		long double distance = (to_variant - from_variant) * 0.000001 * ((long double) recomb_rate) * 4.0L * effective_N;
		// use Li-Stephans pair HMM transitions TODO: correct?
		long double recomb_prob = (1.0L - exp(-distance / (long double) nr_paths) )* (1.0L / (long double) nr_paths);
		long double no_recomb_prob = exp(-distance / (long double) nr_paths) + recomb_prob;
		this->probabilities = {no_recomb_prob*no_recomb_prob, no_recomb_prob*recomb_prob, recomb_prob*recomb_prob};
	}
}

long double CollapsedTransitionProbabilityComputer::compute_transition_prob(DynamicBitset a_id1, DynamicBitset a_id2, DynamicBitset a_id3, DynamicBitset a_id4){
	DynamicBitset hap_1 = a_id1 & a_id3;
	DynamicBitset hap_2 = a_id2 & a_id4;
	size_t count_hap1 = hap_1.count();
	size_t count_hap2 = hap_2.count();
	size_t prod_hap1 = a_id1.count() * a_id3.count() - count_hap1;
	size_t prod_hap2 = a_id2.count() * a_id4.count() - count_hap2;
	return (count_hap1*count_hap2*this->probabilities[0])   +   ( (prod_hap1*count_hap2 + prod_hap2*count_hap1)*this->probabilities[1] )  +  (prod_hap1*prod_hap2*this->probabilities[2]);
}

long double CollapsedTransitionProbabilityComputer::compute_transition_start(DynamicBitset a_id1, DynamicBitset a_id2) {
	return 1.0L * a_id1.count() * a_id2.count();
}
