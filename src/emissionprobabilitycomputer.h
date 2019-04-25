#ifndef EMISSIONPROBABILITYCOMPUTER_H
#define EMISSIONPROBABILITYCOMPUTER_H

#include <vector>
#include <string>
#include <map>
#include "copynumber.h"

class EmissionProbabilityComputer {
public:
	EmissionProbabilityComputer(size_t nr_paths);
	void insert_kmer(std::string kmer, CopyNumber cn, std::vector<int> paths);
	long double get_emission_probability(int path1, int path2);
	bool contains_kmer(std::string kmer);
	CopyNumber get_copynumber_of(std::string kmer);
	std::vector<int> get_paths_of(std::string kmer);

private:
	size_t nr_paths;
	size_t current_index;
	std::map<std::string, size_t> kmer_to_index;
	std::map<size_t,CopyNumber> kmer_to_copynumber;
	std::map<size_t,std::vector<bool> > kmer_to_path;
};
# endif // EMISSIONPROBABILITYCOMPUTER_H
