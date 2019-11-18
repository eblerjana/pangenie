#include <cassert>
#include <math.h>
#include "collapsedtransitionprobabilitycomputer.hpp"
#include <iostream>
#include "dynamicbitset.hpp"

using namespace std;
// size_t from_variant, size_t to_variant, const UniqueKmers* uk_from, const UniqueKmers* uk_to, double recomb_rate, bool uniform = false, long double effective_N = 25000.0L
CollapsedTransitionProbabilityComputer::CollapsedTransitionProbabilityComputer(size_t from_variant, size_t to_variant, const UniqueKmers* uk_from, const UniqueKmers* uk_to, double recomb_rate, size_t nr_paths, bool uniform, long double effective_N)
	:uk_from(uk_from),
	 uk_to(uk_to)
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

long double CollapsedTransitionProbabilityComputer::compute_transition_prob(size_t allele_id1, size_t allele_id2, size_t allele_id3, size_t allele_id4){
	DynamicBitset a_id1 = this->uk_from->get_paths_of_allele(allele_id1);
	DynamicBitset a_id2 = this->uk_from->get_paths_of_allele(allele_id2);
	DynamicBitset a_id3 = this->uk_to->get_paths_of_allele(allele_id3);
	DynamicBitset a_id4 = this->uk_to->get_paths_of_allele(allele_id4);
	DynamicBitset hap_1 = a_id1 & a_id3;
	DynamicBitset hap_2 = a_id2 & a_id4;
	size_t count_hap1 = hap_1.count();
	size_t count_hap2 = hap_2.count();
	size_t prod_hap1 = a_id1.count() * a_id3.count() - count_hap1;
	size_t prod_hap2 = a_id2.count() * a_id4.count() - count_hap2;
	
	return (count_hap1*count_hap2*this->probabilities[0])   +   ( (prod_hap1*count_hap1 + prod_hap2*count_hap2)*this->probabilities[1] )  +  (prod_hap1*prod_hap2*this->probabilities[2]);
}
