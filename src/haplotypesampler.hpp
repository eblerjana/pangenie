#ifndef HAPLOTYPE_SAMPLER_HPP
#define HAPLOTYPE_SAMPLER_HPP

#include <vector>
#include <memory>
#include <stdexcept>
#include <cassert>
#include "uniquekmers.hpp"


struct DPColumn {
	std::vector<unsigned int> column;
};

struct SampledPaths {
	std::vector<std::vector<size_t>> sampled_paths;
	/**
	* Given a column index, and a vector, mark indexes
	* occuring in this column as False.  
	**/
	std::vector<bool> mask_indexes(size_t column_index, size_t max_index) {
		std::vector<bool> masked(max_index+1, true);
		for (auto p : sampled_paths) {
			if (column_index >= p.size()) {
				throw std::runtime_error("HaplotypeSampler::SampledPaths::mask_indexes: column_index exceeds number of columns.");
			}
			size_t index = p[column_index];
			if (index > max_index) {
				throw std::runtime_error("HaplotypeSampler::SampledPaths::mask_indexes: observed index exceeds max_index.");
			}
			masked[index] = false;
		}
		return masked;
	}
};


class HaplotypeSampler {
public:
	/**
	* @param unique_kmers stores the set of unique kmers for each variant position
	* @param size size of the subsampled graph. Viterbi will be run this number of times
	**/
	HaplotypeSampler(std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers, size_t size);
	void rank_haplotypes() const;

	// keeping it public for testing purposes ..
	void get_column_minima(std::vector<unsigned int>& column, std::vector<bool>& mask, size_t& first_id, size_t& second_id, unsigned int& first_val, unsigned int& second_val);

	
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

	/** Update the UniqueKmers object to represent the downsampled graph **/
	void update_unique_kmers();

	std::vector<std::shared_ptr<UniqueKmers>>* unique_kmers;
	std::vector<DPColumn*> viterbi_columns;
	SampledPaths sampled_paths;
	std::vector< std::vector<size_t>* > viterbi_backtrace_columns;
	std::vector<bool> prev_indexes;

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
