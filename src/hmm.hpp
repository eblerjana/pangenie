#ifndef HMM_H
#define HMM_H

#include <vector>
#include "uniquekmers.hpp"
#include "columnindexer.hpp"
#include "transitionprobabilitycomputer.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "probabilitytable.hpp"

/** Respresents the genotyping HMM. **/

class HMM {
public:
	/** 
	* @param unique_kmers stores the set of unique kmers for each variant position.
	* @param run_genotyping run genotyping (Forward backward)
	* @param run_phasing run phasing (Viterbi)
	* @param recombrate recombination rate
	* @param uniform use uniform transition probabilities
	* @param effective_N effective population size
	* @param only_paths only use these paths and ignore others that might be in unique_kmers.
	**/
	HMM(std::vector<UniqueKmers*>* unique_kmers, ProbabilityTable* probabilities, bool run_genotyping, bool run_phasing, double recombrate = 1.26, bool uniform = false, long double effective_N = 25000.0L, std::vector<unsigned short>* only_paths = nullptr, bool normalize = true);
	std::vector<GenotypingResult> get_genotyping_result() const;
	/** moves the GenotypingResults to the caller such that they will no longer be stored in the class. Use with care! **/
	std::vector<GenotypingResult> move_genotyping_result();
	~HMM();

private:
	std::vector<ColumnIndexer*> column_indexers;
	std::vector< std::vector<long double>* > forward_columns;
	std::vector< long double > forward_normalization_sums;
	std::vector<long double>* previous_backward_column;
	std::vector< std::vector<long double>* > viterbi_columns;
	std::vector<UniqueKmers*>* unique_kmers;
	ProbabilityTable* probabilities;
	std::vector< std::vector<size_t>* > viterbi_backtrace_columns;
	std::vector< GenotypingResult > genotyping_result;
	double recombrate;
	bool uniform;
	long double effective_N;
	void compute_forward_prob();
	void compute_backward_prob();
	void compute_viterbi_path();
	void index_columns(std::vector<unsigned short>* only_paths = nullptr);
	void compute_forward_column(size_t column_index);
	void compute_backward_column(size_t column_index);
	void compute_viterbi_column(size_t column_index);

	template<class T>
	void init(std::vector< T* >& c, size_t size) {
		for (size_t i = 0; i < c.size(); ++i) {
			if (c[i] != nullptr) delete c[i];
		}
		c.assign(size, nullptr);
	}
};

#endif // HMM_H
