#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <cassert>
#include "emissionprobabilitycomputer.hpp"

using namespace std;

EmissionProbabilityComputer::EmissionProbabilityComputer(shared_ptr<UniqueKmers> uniquekmers, ProbabilityTable* probabilities)
	:uniquekmers(uniquekmers),
	 probabilities(probabilities),
	 all_zeros(true)
{
	vector<unsigned char> unique_alleles;
	uniquekmers->get_allele_ids(unique_alleles);
	unsigned char max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
	this->state_to_prob = vector< vector<long double>>(max_allele+1, vector<long double>(max_allele+1));

	for (auto a1 : unique_alleles) {
		for (auto a2 : unique_alleles) {
			bool a1_is_undefined = uniquekmers->is_undefined_allele(a1);
			bool a2_is_undefined = uniquekmers->is_undefined_allele(a2);
			this->state_to_prob[a1][a2] = compute_emission_probability(a1, a2, a1_is_undefined, a2_is_undefined);
			if (this->state_to_prob[a1][a2] > 0) this->all_zeros = false;
		}
	}

//	if (this->all_zeros) cerr << "EmissionProbabilities at position " << uniquekmers->get_variant_position() << " are all zero. Set to uniform." << endl;
}

long double EmissionProbabilityComputer::get_emission_probability(unsigned char allele_id1, unsigned char allele_id2) const {
	if (this->all_zeros) return 1.0L;
	return this->state_to_prob[allele_id1][allele_id2];
}

long double EmissionProbabilityComputer::compute_emission_probability(unsigned char allele_id1, unsigned char allele_id2, bool a1_undefined, bool a2_undefined){
	long double result = 1.0L;
	for (size_t i = 0; i < this->uniquekmers->size(); ++i){
		unsigned int expected_kmer_count = this->uniquekmers->alleles.at(allele_id1).first.get_position(i) + this->uniquekmers->alleles.at(allele_id2).first.get_position(i);
		if (a1_undefined && a2_undefined) {
			// all kmers can have copy numbers 0-2
			result *= (1.0L / 3.0L) * (this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(0) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(1) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(2));
		} else if (a1_undefined || a2_undefined) {
			// two possible copy numbers
			assert (expected_kmer_count < 2);
			result *= 0.5L * (this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count) + this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count + 1));
		} else {
			// expected kmer count is known
			result *= this->probabilities->get_probability(this->uniquekmers->local_coverage, this->uniquekmers->kmer_to_count[i]).get_probability_of(expected_kmer_count);
		}
	}
	return result;
}
