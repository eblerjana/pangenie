#include <map>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "samplingemissions.hpp"

using namespace std;

SamplingEmissions::SamplingEmissions(shared_ptr<UniqueKmers> unique_kmers) {
	vector<unsigned char> unique_alleles;
	unique_kmers->get_allele_ids(unique_alleles);
	unsigned char max_allele = *max_element(std::begin(unique_alleles), std::end(unique_alleles));
	this->allele_penalties = vector<unsigned char>(max_allele+1);
	unsigned int default_penalty = 25;

	// precompute penalities
	map<unsigned char, float> fractions = unique_kmers->covered_kmers_on_alleles();
	for (auto it = fractions.begin(); it != fractions.end(); ++it){
		// penalize undefined alleles the same way as uncovered alleles
		if (unique_kmers->is_undefined_allele(it->first)) {
			this->allele_penalties[it->first] = default_penalty;
			continue;
		}

		if (it->second > 0.0) {
			this->allele_penalties[it->first] = -10.0 * log10(it->second);
			assert(this->allele_penalties[it->first] < default_penalty);
		} else {
			// Note: this value is based on the max number of unique kmers (which is 300)
			this->allele_penalties[it->first] = default_penalty;
		}
	}
}

unsigned int SamplingEmissions::get_emission_cost(unsigned char allele_id) {
	return this->allele_penalties[allele_id];
}
