#include <algorithm>
#include <cassert>
#include <math.h>
#include "genotypingresult.hpp"

using namespace std;

GenotypingResult::GenotypingResult()
	:haplotype_1(0),
	 haplotype_2(0),
	 nr_unique_kmers(0)
{}

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
	this->genotype_to_likelihood[genotype] += value;
}

void GenotypingResult::add_first_haplotype_allele(unsigned char allele) {
	this->haplotype_1 = allele;
}

void GenotypingResult::add_second_haplotype_allele(unsigned char allele) {
	this->haplotype_2 = allele;
}

long double GenotypingResult::get_genotype_likelihood (unsigned char allele1, unsigned char allele2) const {
	pair<unsigned char, unsigned char> genotype = genotype_from_alleles(allele1, allele2);
	auto it = this->genotype_to_likelihood.find(genotype);
	if (it != this->genotype_to_likelihood.end()) {
		return this->genotype_to_likelihood.at(genotype);
	} else {
		return 0.0L;
	}
}

vector<long double> GenotypingResult::get_all_likelihoods (size_t nr_alleles) const {
	assert (nr_alleles < 256);
	// determine number of possible genotypes
	size_t nr_genotypes = (nr_alleles * (nr_alleles + 1)) / 2;
	vector<long double> result(nr_genotypes, 0.0L);
	for (auto const& l : this->genotype_to_likelihood) {
		unsigned char allele1 = l.first.first;
		unsigned char allele2 = l.first.second;

		// determine index (according to VCF-specification)
		size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
		if (index >= result.size()) {
			throw runtime_error("GenotypeResult::get_all_likelihoods: genotype does not match number of alleles.");
		}
		result[index] = l.second;
	}
	return result;
}

size_t GenotypingResult::get_genotype_quality (unsigned char allele1, unsigned char allele2) const {
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

pair<unsigned char, unsigned char> GenotypingResult::get_haplotype() const {
	return make_pair(this->haplotype_1, this->haplotype_2);
}

void GenotypingResult::divide_likelihoods_by(long double value) {
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		it->second = it->second / value;
	}
}

pair<int, int> GenotypingResult::get_likeliest_genotype() const {
	long double best_value = 0.0L;
	pair<unsigned char, unsigned char> best_genotype = make_pair(0,0); 
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

	return best_genotype;
}

size_t GenotypingResult::get_nr_unique_kmers() const {
	return this->nr_unique_kmers;
}

size_t GenotypingResult::set_nr_unique_kmers(size_t nr_unique) {
	this->nr_unique_kmers = nr_unique;
}

ostream& operator<<(ostream& os, const GenotypingResult& res) {
	os << "haplotype allele 1: " << res.haplotype_1 << endl;
	os << "haplotype allele 2: " << res.haplotype_2 << endl;
	for (auto it = res.genotype_to_likelihood.begin(); it != res.genotype_to_likelihood.end(); ++it) {
		os << it->first.first << "/" << it->first.second << ": " << it->second << endl;
	}
	return os;
}
