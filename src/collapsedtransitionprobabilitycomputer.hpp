#ifndef COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP
#define COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP

#include "uniquekmers.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class CollapsedTransitionProbabilityComputer {
public:
	CollapsedTransitionProbabilityComputer(size_t from_variant, size_t to_variant, const UniqueKmers* uk_from, const UniqueKmers* uk_to, double recomb_rate, bool uniform = false, long double effective_N = 25000.0L);
	long double compute_transition_prob(size_t allele_id1, size_t allele_id2, size_t allele_id3, size_t allele_id4);
private:
	std::vector<long double> probabilities;
	bool uniform;
	const UniqueKmers* uk_from;
	const UniqueKmers* uk_to;
	
};
#endif // COLLAPSED_TRANSITIONPROBABILITYCOMPUTER_HPP
