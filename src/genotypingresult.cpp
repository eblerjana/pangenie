#include <algorithm>
#include <math.h>
#include "genotypingresult.hpp"

using namespace std;

GenotypingResult::GenotypingResult()
	:haplotype_1(0),
	 haplotype_2(0)
{}

pair<size_t, size_t> genotype_from_alleles (size_t allele1, size_t allele2) {
	// always put allele with smaller index first
	pair<size_t,size_t> genotype;
	if (allele1 < allele2) {
		genotype = make_pair(allele1, allele2);
	} else {
		genotype = make_pair(allele2, allele1);
	}
	return genotype;
}

void GenotypingResult::add_to_likelihood(size_t allele1, size_t allele2, long double value) {
	pair<size_t, size_t> genotype = genotype_from_alleles(allele1, allele2);
	this->genotype_to_likelihood[genotype] += value;
}

void GenotypingResult::add_first_haplotype_allele(size_t allele) {
	this->haplotype_1 = allele;
}

void GenotypingResult::add_second_haplotype_allele(size_t allele) {
	this->haplotype_2 = allele;
}

long double GenotypingResult::get_genotype_likelihood (size_t allele1, size_t allele2) const {
	pair<size_t, size_t> genotype = genotype_from_alleles(allele1, allele2);
	auto it = this->genotype_to_likelihood.find(genotype);
	if (it != this->genotype_to_likelihood.end()) {
		return this->genotype_to_likelihood.at(genotype);
	} else {
		return 0.0L;
	}
}

vector<long double> GenotypingResult::get_all_likelihoods (size_t nr_alleles) const {
	// determine number of possible genotypes
	size_t nr_genotypes = (nr_alleles * (nr_alleles + 1)) / 2;
	vector<long double> result(nr_genotypes, 0.0L);
	for (auto const& l : this->genotype_to_likelihood) {
		size_t allele1 = l.first.first;
		size_t allele2 = l.first.second;

		// determine index (according to VCF-specification)
		size_t index = ((allele2 * (allele2 + 1)) / 2) + allele1;
		if (index >= result.size()) {
			throw runtime_error("GenotypeResult::get_all_likelihoods: genotype does not match number of alleles.");
		}
		result[index] = l.second;
	}
	return result;
}

size_t GenotypingResult::get_genotype_quality (size_t allele1, size_t allele2) const {
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

pair<size_t, size_t> GenotypingResult::get_haplotype() const {
	return make_pair(this->haplotype_1, this->haplotype_2);
}

void GenotypingResult::divide_likelihoods_by(long double value) {
	for (auto it = this->genotype_to_likelihood.begin(); it != this->genotype_to_likelihood.end(); ++it) {
		it->second = it->second / value;
	}
}

pair<int, int> GenotypingResult::get_likeliest_genotype() const {
	long double best_value = 0.0L;
	pair<size_t,size_t> best_genotype = make_pair(0,0); 
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

ostream& operator<<(ostream& os, const GenotypingResult& res) {
	os << "haplotype allele 1: " << res.haplotype_1 << endl;
	os << "haplotype allele 2: " << res.haplotype_2 << endl;
	for (auto it = res.genotype_to_likelihood.begin(); it != res.genotype_to_likelihood.end(); ++it) {
		os << it->first.first << "/" << it->first.second << ": " << it->second << endl;
	}
	return os;
}
