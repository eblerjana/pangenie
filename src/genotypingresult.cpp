#include <algorithm>
#include <cassert>
#include <math.h>
#include "genotypingresult.hpp"

using namespace std;

GenotypingResult::GenotypingResult(size_t nr_alleles)
	:haplotype_1(0),
	 haplotype_2(0),
	 local_coverage(0),
	 unique_kmers(0)
{
	if (nr_alleles < 2) {
		this->nr_alleles = 2;
	} else {
		this->nr_alleles = nr_alleles;
	}
	
	this->genotype_to_likelihood = vector<long double>((this->nr_alleles * (this->nr_alleles + 1)) / 2, 0.0L);
}

pair<unsigned char, unsigned char> genotype_from_alleles (unsigned char allele1, unsigned char allele2) {
	// always put allele with smaller index first
	pair<unsigned char,unsigned char> genotype;
	if (allele1 < allele2) {
		genotype = make_pair(allele1, allele2);
	} else {
		genotype = make_pair(allele2, allele1);
	}
	return genotype;
}

void GenotypingResult::add_to_likelihood(unsigned char allele1, unsigned char allele2, long double value) {
	pair<unsigned char, unsigned char> genotype = genotype_from_alleles(allele1, allele2);
	// determine index
	size_t index = ((genotype.second * (genotype.second + 1)) / 2) + genotype.first;
	this->genotype_to_likelihood[index] += value;
}

void GenotypingResult::add_first_haplotype_allele(unsigned char allele) {
	this->haplotype_1 = allele;
}

void GenotypingResult::add_second_haplotype_allele(unsigned char allele) {
	this->haplotype_2 = allele;
}

long double GenotypingResult::get_genotype_likelihood (unsigned char allele1, unsigned char allele2) const {
	pair<unsigned char, unsigned char> genotype = genotype_from_alleles(allele1, allele2);
	size_t index = ((genotype.second * (genotype.second + 1)) / 2) + genotype.first;
	return this->genotype_to_likelihood.at(index);
}

vector<long double> GenotypingResult::get_all_likelihoods () const {
	assert (nr_alleles < 256);
	return this->genotype_to_likelihood;
}

GenotypingResult GenotypingResult::get_specific_likelihoods (vector<unsigned char>& alleles) const {
	size_t nr_alleles = alleles.size();
	unsigned char max_allele = 0;
	if (alleles.size() > 0) {
		max_allele = *max_element(std::begin(alleles), std::end(alleles));
	}
	GenotypingResult result(max_allele+1);
	assert (nr_alleles < 256);
	long double sum = 0.0L;
	for (unsigned char i = 0; i < nr_alleles; ++i) {
		for (unsigned char j = i; j < nr_alleles; ++j) {
			long double likelihood = this->get_genotype_likelihood(alleles[i], alleles[j]);
			if (this->haplotype_1 == alleles[i]) result.haplotype_1 = i;
			if (this->haplotype_2 == alleles[j]) result.haplotype_2 = j;
			result.add_to_likelihood(i, j, likelihood);
			sum += likelihood;
		}
	}
	if (sum > 0) result.divide_likelihoods_by(sum);
	return result;
}


size_t GenotypingResult::get_genotype_quality (unsigned char allele1, unsigned char allele2) const {
	// check if likelihoods are normalized
	long double sum = 0.0;
	for (const auto& l : this->genotype_to_likelihood) {
		sum += l;
	}

	if (abs(sum-1) > 0.0000000001) {
		throw runtime_error("GenotypingResult::get_genotype_quality: genotype quality can only be computed from normalized likelihoods.");
	}
	
	// compute genotype quality
	long double prob_wrong = 1.0L - get_genotype_likelihood(allele1, allele2);
	if (prob_wrong > 0.0) {
		return (size_t) (-10 * log10(prob_wrong));
	} else {
		// TODO: default value ok?
		return 10000;
	}
}

pair<unsigned char, unsigned char> GenotypingResult::get_haplotype() const {
	return make_pair(this->haplotype_1, this->haplotype_2);
}


void GenotypingResult::divide_likelihoods_by(long double value) {
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		*it = (*it / value);
	}
}


pair<int, int> GenotypingResult::get_likeliest_genotype() const {

	long double best_value = 0.0L;
	pair<unsigned char, unsigned char> best_genotype = make_pair(0,0);

	for (unsigned char allele1 = 0; allele1 < this->nr_alleles; ++allele1) {
		for (unsigned char allele2 = allele1; allele2 < this->nr_alleles; ++allele2) {
			size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
			if (this->genotype_to_likelihood.at(index) >= best_value) {
				best_value = this->genotype_to_likelihood.at(index);
				best_genotype = make_pair(allele1, allele2);
			}
		}
	}


	for (unsigned char allele1 = 0; allele1 < this->nr_alleles; ++allele1) {
		for (unsigned char allele2 = allele1; allele2 < this->nr_alleles; ++allele2) {
			size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
			if (best_genotype != pair<unsigned char, unsigned char>(allele1, allele2)) {
				if (abs(this->genotype_to_likelihood.at(index) - best_value) < 0.0000000001) {
					// not unique
					return make_pair(-1,-1);
				}
			}
		}
	}

	// if best genotype has likelihood 0 (this can happen if there is only one entry), return ./.
	if (best_value > 0.0L) {
		return best_genotype;
	} else {
		return make_pair(-1,-1);
	}
}

ostream& operator<<(ostream& os, const GenotypingResult& res) {
	os << "haplotype allele 1: " << res.haplotype_1 << endl;
	os << "haplotype allele 2: " << res.haplotype_2 << endl;
	os << "local coverage: " << res.local_coverage << endl;
	os << "nr of unique kmers: " << res.unique_kmers << endl;

	for (size_t allele1 = 0; allele1 < res.nr_alleles; ++allele1) {
		for (size_t allele2 = allele1; allele2 < res.nr_alleles; ++allele2) {
			size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
			os << (unsigned int) allele1 << "/" << allele2 << ":" << res.genotype_to_likelihood.at(index) << endl;
		}
	}
	return os;
}


void GenotypingResult::combine(GenotypingResult& likelihoods) {

	if (this->nr_alleles > likelihoods.nr_alleles) {
		std::transform (this->genotype_to_likelihood.begin(), this->genotype_to_likelihood.begin() + likelihoods.genotype_to_likelihood.size(), likelihoods.genotype_to_likelihood.begin(), this->genotype_to_likelihood.begin(), std::plus<long double>());
	} else {
		this->nr_alleles = likelihoods.nr_alleles;
		this->genotype_to_likelihood.resize(likelihoods.genotype_to_likelihood.size());
		std::transform (this->genotype_to_likelihood.begin(), this->genotype_to_likelihood.end(), likelihoods.genotype_to_likelihood.begin(), this->genotype_to_likelihood.begin(), std::plus<long double>());
	}
}

void GenotypingResult::normalize () {
	// sum up probabilities
	long double normalization_sum = 0.0L;
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		normalization_sum += *it;
	}

	if (normalization_sum > 0) {
		this->divide_likelihoods_by(normalization_sum);
	}
}


void GenotypingResult::set_unique_kmers(unsigned short nr_unique_kmers) {
	assert (nr_unique_kmers < 65536);
	this->unique_kmers = nr_unique_kmers;
}


void GenotypingResult::set_coverage(unsigned short coverage) {
	assert (coverage < 65536);
	this->local_coverage = coverage;
}


unsigned short GenotypingResult::nr_unique_kmers() const {
	return this->unique_kmers;
}


unsigned short GenotypingResult::coverage() const {
	return this->local_coverage;
}

size_t GenotypingResult::size() const {
	return this->nr_alleles;
}
