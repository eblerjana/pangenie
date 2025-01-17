#include <algorithm>
#include <cassert>
#include <math.h>
#include "genotypingresult.hpp"

using namespace std;

GenotypingResult::GenotypingResult()
	:haplotype_1(0),
	 haplotype_2(0),
	 local_coverage(0),
	 unique_kmers(0)
{}

pair<unsigned short, unsigned short> genotype_from_alleles (unsigned short allele1, unsigned short allele2) {
	// always put allele with smaller index first
	pair<unsigned short,unsigned short> genotype;
	if (allele1 < allele2) {
		genotype = make_pair(allele1, allele2);
	} else {
		genotype = make_pair(allele2, allele1);
	}
	return genotype;
}

void GenotypingResult::add_to_likelihood(unsigned short allele1, unsigned short allele2, long double value) {
	pair<unsigned short, unsigned short> genotype = genotype_from_alleles(allele1, allele2);
	this->genotype_to_likelihood[genotype] += value;
}

void GenotypingResult::add_first_haplotype_allele(unsigned short allele) {
	this->haplotype_1 = allele;
}

void GenotypingResult::add_second_haplotype_allele(unsigned short allele) {
	this->haplotype_2 = allele;
}

long double GenotypingResult::get_genotype_likelihood (unsigned short allele1, unsigned short allele2) const {
	pair<unsigned short, unsigned short> genotype = genotype_from_alleles(allele1, allele2);
	auto it = this->genotype_to_likelihood.find(genotype);
	if (it != this->genotype_to_likelihood.end()) {
		return this->genotype_to_likelihood.at(genotype);
	} else {
		return 0.0L;
	}
}

vector<long double> GenotypingResult::get_all_likelihoods (size_t nr_alleles) const {
	assert (nr_alleles < 65536);

	// determine number of possible genotypes
	size_t nr_genotypes = (nr_alleles * (nr_alleles + 1)) / 2;

	vector<long double> result(nr_genotypes, 0.0L);
	for (auto const& l : this->genotype_to_likelihood) {
		unsigned short allele1 = l.first.first;
		unsigned short allele2 = l.first.second;

		// determine index (according to VCF-specification)
		size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
		if (index >= result.size()) {
			throw runtime_error("GenotypeResult::get_all_likelihoods: genotype does not match number of alleles.");
		}
		result[index] = l.second;
	}
	return result;
}

GenotypingResult GenotypingResult::get_specific_likelihoods (vector<unsigned short>& alleles) const {
	GenotypingResult result;
	size_t nr_alleles = alleles.size();
	assert (nr_alleles < 65536);
	long double sum = 0.0L;
	for (unsigned short i = 0; i < nr_alleles; ++i) {
		for (unsigned short j = i; j < nr_alleles; ++j) {
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

size_t GenotypingResult::get_genotype_quality (unsigned short allele1, unsigned short allele2) const {
	// check if likelihoods are normalized
	long double sum = 0.0;
	for (const auto& l : this->genotype_to_likelihood) {
		sum += l.second;
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

pair<unsigned short, unsigned short> GenotypingResult::get_haplotype() const {
	return make_pair(this->haplotype_1, this->haplotype_2);
}

void GenotypingResult::divide_likelihoods_by(long double value) {
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		it->second = it->second / value;
	}
}

pair<int, int> GenotypingResult::get_likeliest_genotype() const {
	// if empty, set genotype to unknown
	if (this->genotype_to_likelihood.size() == 0) {
		return make_pair(-1,-1);
	}

	long double best_value = 0.0L;
	pair<unsigned short, unsigned short> best_genotype = make_pair(0,0); 
	for (auto const& l : this->genotype_to_likelihood) {
		if (l.second >= best_value) {
			best_value = l.second;
			best_genotype = l.first;
		}
	}

	// make sure there is a unique maximum
	for (auto const& l : this->genotype_to_likelihood) {
		if (best_genotype != l.first) {
			if (abs(l.second-best_value) < 0.0000000001) {
				// not unique
				return make_pair(-1,-1);
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
	for (auto it = res.genotype_to_likelihood.begin(); it != res.genotype_to_likelihood.end(); ++it) {
		os << (unsigned int) it->first.first << "/" << (unsigned int) it->first.second << ": " << it->second << endl;
	}
	return os;
}

void GenotypingResult::combine(GenotypingResult& likelihoods) {
	for (auto it = likelihoods.genotype_to_likelihood.begin(); it != likelihoods.genotype_to_likelihood.end(); ++it) {
		pair<unsigned short, unsigned short> genotype = it->first;
		this->genotype_to_likelihood[genotype] += it->second;
	}
}

void GenotypingResult::normalize () {
	// sum up probabilities
	long double normalization_sum = 0.0L;
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		normalization_sum += it->second;
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

bool GenotypingResult::contains_no_likelihoods() const {
	return this->genotype_to_likelihood.empty();
}