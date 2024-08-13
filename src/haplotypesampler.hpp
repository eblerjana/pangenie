#ifndef HAPLOTYPE_SAMPLER_HPP
#define HAPLOTYPE_SAMPLER_HPP

#include <vector>
#include <memory>
#include "uniquekmers.hpp"


struct DPColumn {
	std::vector<long double> column;
};

struct SampledPaths {
	std::vector<std::vector<size_t>> sampled_paths;
};


class HaplotypeSampler {
public:
	/**
	* @param unique_kmers stores the set of unique kmers for each variant position
	* @param size size of the subsampled graph. Viterbi will be run this number of times
	**/
	HaplotypeSampler(std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers, size_t size);
	void rank_haplotypes() const;

	
private:
	/** Do one Viterbi pass and store the paths that have been used. 
	* This function can be applied several times and keeps track of used Viterbi path,
	* so that these nodes are ignored in the next pass (= call to this function). Used Viterbi
	* paths are stored in member variable "sampled_paths".
	**/
	void compute_viterbi_path();
	
	/**
	* Compute one column. Make sure that when iterating through the paths (= states), those paths
	* that were already used in previous Viterbi runs, are ignored. For this purpose, determine a list
	* of respective path ids first, so that these can be ignored while iterating the states.
	**/
	void compute_viterbi_column(size_t column_index);

	/** Update the UniqueKmers object to represent the downsampled grah **/
	void update_unique_kmers();

	std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers;
	std::vector<DPColumn*> viterbi_columns;
	SampledPaths sampled_paths;
	std::vector< std::vector<size_t>* > viterbi_backtrace_columns;

	template<class T>
	void init(std::vector< T* >& c, size_t size) {
		for (size_t i = 0; i < c.size(); ++i) {
			if (c[i] != nullptr) delete c[i];
		}
		c.assign(size, nullptr);
	}


};

#endif // HAPLOTYPE_SAMPLER_HPP




/**
* Implementation ideas:
* - in each iteration of Viterbi, maxima need to be computed
* - this can be done efficiently by precomputing maxima for all sets: {p1, ..., pi-1, pi+1, ..., pn} (excluding each pi): 
*   compute overall maximum: x = max(p1,...,pn). Compute y = max(p1,...,px-1,px+1,...,pn). precomputed max for all sets with px is x, for all others y    
* - for each cell i compute: max(t0 * pi, t1*max(p1,...,pn)), the max(p1,..pn) term is precomputed
*
* Need: one class/function computing the transition probability / score /penalty, and one for the Emission.
* Transitions: compute phred-scaled Li stephans / WH transition score
* Emissions / local cost: 
**/
