#ifndef COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP
#define COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP

#include "dynamicbitset.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class CollapsedTransitionProbabilityComputer {
public:
	CollapsedTransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, size_t nr_paths, bool uniform = false, long double effective_N = 25000.0L);
	long double compute_transition_prob(DynamicBitset paths_allele1, DynamicBitset paths_allele2, DynamicBitset paths_allele3, DynamicBitset paths_allele4);
	long double compute_transition_start(DynamicBitset paths_allele1, DynamicBitset paths_allele2);
private:
	std::vector<long double> probabilities;
	size_t nr_paths;
	
};
#endif // COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP
