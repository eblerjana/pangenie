#ifndef PROBABILITYTABLE_HPP
#define PROBABILITYTABLE_HPP

#include <vector>
#include "copynumber.hpp"
#include <iostream>

/** 
* Pre-computes probabilities for kmer copy numbers and read kmer counts.
**/

class ProbabilityTable {
public:
	ProbabilityTable();
	ProbabilityTable(unsigned short cov_min, unsigned short cov_max, unsigned short count_max, long double regularization_const);
	CopyNumber get_probability (unsigned short kmer_coverage, unsigned short read_kmer_count) const;
	/** function can be used to modify probabilities stored in the table. Mainly used for testing purposes. **/
	void modify_probability(unsigned short kmer_coverage, unsigned short read_kmer_count, CopyNumber prob);
	friend std::ostream& operator<<(std::ostream& os, const ProbabilityTable& table);
private:
	unsigned short cov_min;
	unsigned short cov_max;
	unsigned short count_max;
	long double regularization_const;
	std::vector<std::vector<CopyNumber>> probabilities;
	long double poisson(long double mean, unsigned int value) const;
	long double geometric(long double p, unsigned int value) const;
	CopyNumber compute_probability (unsigned short kmer_coverage, unsigned short read_kmer_count) const;
};
#endif // PROBABILITYTABLE_HPP
