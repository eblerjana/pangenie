#ifndef PATH_SAMPLER_HPP
#define PATH_SAMPLER_HPP

#include <vector>

class PathSampler {
public:
	PathSampler(unsigned short total_number);
	// select one random subset of paths
	void select_single_subset(std::vector<unsigned short>& result, unsigned short sample_size) const;
	// select random subset of paths n times
	void select_multiple_subsets(std::vector<std::vector<unsigned short>>& result, unsigned short sample_size, unsigned short n) const;
	// partition set of paths
	void partition_paths(std::vector<std::vector<unsigned short>>& result, unsigned short sample_size) const;
	// partition paths in sample aware manner
	void partition_samples(std::vector<std::vector<unsigned short>>& result, unsigned short sample_size) const;
private:
	unsigned short total_number; 
};


#endif // PATH_SAMPLER_HPP

