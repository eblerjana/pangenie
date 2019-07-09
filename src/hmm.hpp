#ifndef HMM_H
#define HMM_H

#include <vector>
#include "uniquekmers.hpp"
#include "columnindexer.hpp"
#include "transitionprobabilitycomputer.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"

/** Respresents the genotyping HMM. **/

class HMM {
public:
	/** 
	* @param unique_kmers stores the set of unique kmers for each variant position.
	* @param variants list of all variants to consider
	* @param recombrate recombination rate
	**/
	HMM(std::vector<UniqueKmers>* unique_kmers, const std::vector<Variant>& variants, double recombrate = 1.26);
	const std::vector<GenotypingResult>& get_genotyping_result() const;
	~HMM();

private:
	std::vector<ColumnIndexer*> column_indexers;
	std::vector< std::vector<long double>* > forward_columns;
	std::vector<long double>* previous_backward_column;
	std::vector< std::vector<long double>* > viterbi_columns;
	std::vector<UniqueKmers>* unique_kmers;
	std::vector<TransitionProbabilityComputer*> transition_prob_computers;
	std::vector< std::vector<size_t>* > viterbi_backtrace_columns;
	std::vector< GenotypingResult > genotyping_result;
	void compute_forward_prob();
	void compute_backward_prob(const std::vector<Variant>& variants);
	void compute_viterbi_path(const std::vector<Variant>& variants);
	void index_columns();
	void compute_forward_column(size_t column_index);
	void compute_backward_column(size_t column_index, const std::vector<Variant>& variants);
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