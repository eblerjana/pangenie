#ifndef HAPLOTYPE_SAMPLER_HPP
#define HAPLOTYPE_SAMPLER_HPP

#include <vector>
#include <memory>
#include "uniquekmers.hpp"


class HaplotypeSampler {
public:
	HaplotypeSampler(std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers);
	void rank_haplotypes() const;
private:
	std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers;
};

#endif // HAPLOTYPE_SAMPLER_HPP