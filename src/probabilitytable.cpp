#include "probabilitytable.hpp"
#include <math.h>
#include <stdexcept>

using namespace std;

double get_error_param(double kmer_coverage) {
	double cn0;
	if (kmer_coverage < 10.0) {
		cn0 = 0.99;
	} else if (kmer_coverage < 20) {
		cn0 = 0.95;
	} else if (kmer_coverage < 40) {
		cn0 = 0.9;
	} else {
		cn0 = 0.8;
	}
	return cn0;
}

ProbabilityTable::ProbabilityTable()
	:cov_min(0),
	 cov_max(0),
	 count_max(0),
	 regularization_const(0.0L)
{}

ProbabilityTable::ProbabilityTable(unsigned short cov_min, unsigned short cov_max, unsigned short count_max, long double regularization_const)
	:cov_min(cov_min),
	 cov_max(cov_max),
	 count_max(count_max),
	 regularization_const(regularization_const)
{
	// initialize table
	for (unsigned short i = 0; i < count_max; ++i) {
		this->probabilities.push_back(vector<CopyNumber>(cov_max-cov_min));
	}

	for (unsigned short i = 0; i < this->count_max; ++i) {
		// precompute probabilities for each read kmer count
		for (unsigned short j = 0; (j + this->cov_min) < this->cov_max; ++j) {
			this->probabilities[i][j] = compute_probability(j + this->cov_min, i);
		}
	}
}

CopyNumber ProbabilityTable::get_probability (unsigned short kmer_coverage, unsigned short read_kmer_count) const {
	if ((kmer_coverage >= this->cov_min) && (kmer_coverage < this->cov_max) && (read_kmer_count < this->count_max)) {
		return this->probabilities.at(read_kmer_count).at(kmer_coverage - this->cov_min);
	} else {
		return compute_probability(kmer_coverage, read_kmer_count);
	}
}

CopyNumber ProbabilityTable::compute_probability(unsigned short kmer_coverage, unsigned short read_kmer_count) const {
			long double p_cn0 = geometric(get_error_param(kmer_coverage), read_kmer_count);
			long double p_cn1 = poisson(kmer_coverage / 2.0, read_kmer_count);
			long double p_cn2 = poisson(kmer_coverage, read_kmer_count);

			if (this->regularization_const > 0) {
				return CopyNumber(p_cn0, p_cn1, p_cn2, regularization_const);
			} else {
				return CopyNumber(p_cn0, p_cn1, p_cn2);
			}
}

void ProbabilityTable::modify_probability(unsigned short kmer_coverage, unsigned short read_kmer_count, CopyNumber prob) {
	if ((kmer_coverage >= this->cov_min) && (kmer_coverage < this->cov_max) && (read_kmer_count < this->count_max)) {
		this->probabilities[read_kmer_count][kmer_coverage - this->cov_min] = prob;
	} else {
		throw runtime_error("ProbabilityTable::modify_probability: no precomputed values for these parameters.");
	}
}

long double ProbabilityTable::poisson(long double mean, unsigned int value) const {
	long double sum = 0.0L;
	int v = (int) value;
	for (size_t i = 1; i <= value; ++i) sum += log(i);
	long double log_val = -mean + v * log(mean) - sum;
	return exp(log_val);
}

long double ProbabilityTable::geometric(long double p, unsigned int value) const {
	return pow(1.0L - p, value)*p;
}

ostream& operator<<(ostream& os, const ProbabilityTable& var) {
	os << "\t";
	for (unsigned short i = var.cov_min; i < var.cov_max; ++i) {
		if (i > var.cov_min) os << "\t\t\t";
		os << i;
	}
	os << "\n";
	for (unsigned short i = 0; i < var.count_max; ++i) {
		os << i << "\t";
		for (unsigned short j = 0; (j + var.cov_min) < var.cov_max; ++j) {
			if (j > 0) os << "\t";
			os << var.probabilities.at(i).at(j).get_probability_of(0) << "\t";
			os << var.probabilities.at(i).at(j).get_probability_of(1) << "\t";
			os << var.probabilities.at(i).at(j).get_probability_of(2);
		}
		os << "\n";
	}
	return os;
}
