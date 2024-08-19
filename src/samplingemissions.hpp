#ifndef SAMPLING_EMISSIONS
#define SAMPLING_EMISSIONS

#include <memory>
#include "uniquekmers.hpp"

class SamplingEmissions {
public:
	SamplingEmissions(std::shared_ptr<UniqueKmers> uniquekmers);
	float get_emission_cost(unsigned char allele_id);
private:
	std::vector<unsigned int> allele_penalties;

};

#endif // SAMPLING_EMISSIONS
