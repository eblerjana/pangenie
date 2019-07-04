#ifndef TRANSITIONPROBABILITYCOMPUTER_HPP
#define TRANSITIONPROBABILITYCOMPUTER_HPP

#include "uniquekmers.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class TransitionProbabilityComputer {
public:
	TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate);
	long double compute_transition_prob(size_t path_id1, size_t path_id2, size_t path_id3, size_t path_id4);
private:
	std::vector<long double> probabilities;
	
};
#endif // TRANSITIONPROBABILITYCOMPUTER_HPP
