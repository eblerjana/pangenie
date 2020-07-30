#ifndef TRANSITIONPROBABILITYCOMPUTER_HPP
#define TRANSITIONPROBABILITYCOMPUTER_HPP

#include "uniquekmers.hpp"

/** 
* Computes the transition probabilities between variants.
**/

class TransitionProbabilityComputer {
public:
	TransitionProbabilityComputer(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, bool uniform = false, long double effective_N = 25000.0L);
	long double compute_transition_prob(unsigned short path_id1, unsigned short path_id2, unsigned short path_id3, unsigned short path_id4);
private:
	std::vector<long double> probabilities;
	bool uniform;
	
};
#endif // TRANSITIONPROBABILITYCOMPUTER_HPP
