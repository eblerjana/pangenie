#ifndef SAMPLING_EMISSIONS
#define SAMPLING_EMISSIONS

#include <memory>
#include "uniquekmers.hpp"

class SamplingEmissions {
public:
	SamplingEmissions(std::shared_ptr<UniqueKmers> uniquekmers);
	unsigned int get_emission_cost(unsigned short allele_id) const;
	void penalize(unsigned short allele_id);
private:
	std::vector<unsigned short> allele_penalties;
	unsigned int default_penalty;

};

#endif // SAMPLING_EMISSIONS
