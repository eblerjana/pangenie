#ifndef PATH_SAMPLER_HPP
#define PATH_SAMPLER_HPP

#include <vector>

class PathSampler {
public:
	PathSampler(size_t total_number);
	// select one random subset of paths
	void select_single_subset(std::vector<size_t>& result, size_t sample_size) const;
	// select random subset of paths n times
	void select_multiple_subsets(std::vector<std::vector<size_t>>& result, size_t sample_size, size_t n) const;
	// partition set of paths
	void partition_paths(std::vector<std::vector<size_t>>& result, size_t sample_size) const;
private:
	size_t total_number; 
};


#endif // PATH_SAMPLER_HPP

