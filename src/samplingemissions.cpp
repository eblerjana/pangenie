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
	this->default_penalty = 25;
	unsigned int undefined_penalty = 50;

	for (auto a : unique_alleles) {
		if (unique_kmers->is_undefined_allele(a)) {
			this->allele_penalties[a] = undefined_penalty;
			continue;
		}

		float fraction = unique_kmers->fraction_present_kmers_on_allele(a);

		if (fraction > 0.0) {
			this->allele_penalties[a] = -10.0 * log10(fraction);
			assert(this->allele_penalties[a] < this->default_penalty);
		} else {
			// Note: this value is based on the max number of unique kmers (which is 300)
			this->allele_penalties[a] = this->default_penalty;
		}
	}
}

unsigned int SamplingEmissions::get_emission_cost(unsigned char allele_id) const {
	return this->allele_penalties[allele_id];
}

void SamplingEmissions::penalize(unsigned char allele_id) {
	this->allele_penalties[allele_id] += 10;
	if (this->allele_penalties[allele_id] > this->default_penalty) {
		// make sure max penality value is at most default + 10
		this->allele_penalties[allele_id] = this->default_penalty + 10;
	}
}