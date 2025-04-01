#ifndef SAMPLING_TRANSITIONS
#define SAMPLING_TRANSITIONS

#include <vector>
#include <cassert>

class SamplingTransitions {

public:
	/**
	* @param from_variant, to_variant positions of variants to compute recominbation cost for
	* @param recomb_rate recombination rate
	* @param number of paths in the panel
	* @param effective_N effective population size
	**/
	SamplingTransitions(size_t from_variant, size_t to_variant, double recomb_rate, unsigned short nr_paths, long double effective_N = 25000.0L);
	/** computes transition cost. If there is a recombination event, cost is defined as the phred-scaled recombination probability
	* if no recombination, cost is zero.
	* @param indicates if recombination event happened between columns
	***/
	unsigned int compute_transition_cost(bool recombination);
private:
	unsigned int cost;
};

#endif // SAMPLING_TRANSITIONS
