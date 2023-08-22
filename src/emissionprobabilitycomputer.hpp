#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "uniquekmers.hpp"
#include "copynumber.hpp"
#include "columnindexer.hpp"
#include "probabilitytable.hpp"

/** 
* Computes the emission probabilities for a variant position.
**/

typedef std::vector<std::vector<long double>> ProbabilityMatrix;

class EmissionProbabilityComputer {
public:
	/**
	* @param uniquekmers all unique kmers for this position
	 **/
	EmissionProbabilityComputer(std::shared_ptr<UniqueKmers> uniquekmers, ProbabilityTable* probabilities);
	/** get emission probability for a state in the HMM **/
	long double get_emission_probability(unsigned char allele_id1, unsigned char allele_id2) const;

private:
	std::shared_ptr<UniqueKmers> uniquekmers;
	ProbabilityTable* probabilities;
	bool all_zeros;
	ProbabilityMatrix state_to_prob;
	long double compute_emission_probability(unsigned char allele1, unsigned char allele2, bool allele1_undefined, bool allele2_undefined);
};
# endif // EMISSIONPROBABILITYCOMPUTER_H
