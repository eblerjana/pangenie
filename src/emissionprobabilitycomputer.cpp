#include <stdexcept>
#include <sstream>
#include <algorithm>
#include "emissionprobabilitycomputer.hpp"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(UniqueKmers* uniquekmers)
	:uniquekmers(uniquekmers),
	 all_zeros(true)
{
	vector<unsigned char> unique_alleles;
	uniquekmers->get_allele_ids(unique_alleles);
	unsigned char max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
	this->state_to_prob = vector< vector<long double>>(max_allele+1, vector<long double>(max_allele+1));

	for (auto a1 : unique_alleles) {
		for (auto a2 : unique_alleles) {
			this->state_to_prob[a1][a2] = compute_emission_probability(a1, a2);
			if (this->state_to_prob[a1][a2] > 0) this->all_zeros = false;
		}
	}

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

long double EmissionProbabilityComputer::get_emission_probability(unsigned char allele_id1, unsigned char allele_id2) const {
	if (this->all_zeros) return 1.0L;
	return this->state_to_prob[allele_id1][allele_id2];
}

long double EmissionProbabilityComputer::compute_emission_probability(unsigned char allele_id1, unsigned char allele_id2){
	long double result = 1.0L;
	// combine the two paths to get expected kmer copy numbers
	CopyNumberAssignment cna = this->uniquekmers->combine_paths(allele_id1, allele_id2);
	for (size_t i = 0; i < this->uniquekmers->size(); ++i){
		unsigned int expected_kmer_count = cna.get_position(i);
		// multiply result by probability of expected kmer count
		result *= uniquekmers->kmer_to_copynumber[i].get_probability_of(expected_kmer_count);
	}
	return result;
}
