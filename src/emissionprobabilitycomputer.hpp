#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include <vector>
#include <string>
#include <unordered_map>
#include "uniquekmers.hpp"
#include "copynumber.hpp"

/** 
* Computes the emission probabilities for a variant position.
**/


struct pair_hash {
	size_t operator() (const std::pair<size_t,size_t> &p) const {
		return ((p.first + p.second) * (p.first + p.second + 1) / 2 + p.second);
	}
};


class EmissionProbabilityComputer {
public:
	/**
	* @param uniquekmers all unique kmers for this position
	 **/
	EmissionProbabilityComputer(UniqueKmers* uniquekmers);
	/** get emission probability for a pair of paths (state in HMM) **/
	long double get_emission_probability(int path1, int path2) const;

private:
	UniqueKmers* uniquekmers;
	std::unordered_map< std::pair<size_t,size_t>, long double, pair_hash> paths_to_prob;
	long double compute_emission_probability(int path1, int path2);
};
# endif // EMISSIONPROBABILITYCOMPUTER_H
