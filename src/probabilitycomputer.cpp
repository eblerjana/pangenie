#include "probabilitycomputer.hpp"
#include <math.h>
#include <stdexcept>

using namespace std;

ProbabilityComputer::ProbabilityComputer() {}

ProbabilityComputer::ProbabilityComputer(long double mean_cn0, long double mean_cn1, long double mean_cn2)
	:means({mean_cn0, mean_cn1, mean_cn2})
{}

void ProbabilityComputer::set_parameters (long double mean_cn0, long double mean_cn1, long  double mean_cn2) {
	if ( (mean_cn0 > 0) && (mean_cn1 > 0) && (mean_cn2 > 0) ) { 
		this->means = {mean_cn0, mean_cn1, mean_cn2};
	} else {
		throw runtime_error("ProbabilityComputer::set_parameters: all means need to be > 0.");
	}
}

long double ProbabilityComputer::get_probability(size_t cn, unsigned int value) const {
	if (this->means.size() == 0) {
		throw runtime_error("ProbabilityComputer::get_probability: parameters have not yet been set. Use set_parameters().");
	}
	if (cn > 2) {
		throw runtime_error("ProbabilityComputer::get_probability: only copynumbers 0, 1 and 2 allowed.");
	}
	if (cn > 0) {
		return poisson(means[cn], value);
	} else {
		return geometric(means[cn], value);
	}
}

long double ProbabilityComputer::poisson(long double mean, unsigned int value) const {
	long double sum = 0.0L;
	int v = (int) value;
	for (size_t i = 1; i <= value; ++i) sum += log(i);
	long double log_val = -mean + v * log(mean) - sum;
	return exp(log_val);
}

long double ProbabilityComputer::geometric(long double p, unsigned int value) const {
	return pow(1.0L - p, value)*p;
}
