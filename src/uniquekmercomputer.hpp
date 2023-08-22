#ifndef UNIQUEKMERCOMPUTER_HPP
#define UNIQUEKMERCOMPUTER_HPP

#include <vector>
#include <string>
#include <memory>
#include "kmercounter.hpp"
#include "graph.hpp"
#include "uniquekmers.hpp"
#include "probabilitytable.hpp"

class UniqueKmerComputer {
public:
	/** 
	* @param genomic_kmers genomic kmer counts
	* @param read_kmers read kmer counts
	* @param variants 
	* @param kmer_coverage needed to compute kmer copy number probabilities
	**/
	UniqueKmerComputer (KmerCounter* genomic_kmers, std::shared_ptr<KmerCounter> read_kmers, std::shared_ptr<Graph> variants, size_t kmer_coverage);
	/** generates UniqueKmers object for each position, ownership of vector is transferred to the caller.
	* @param result	UniqueKmer objects will be stored here
	* @param probabilities pre-computed ProbabilityTable
	* @param delete_processed_variants if set to true, the VariantReader will be motified. Processed Variant objects will be deleted. Use with care!
	 **/
	void compute_unique_kmers(std::vector<std::shared_ptr<UniqueKmers>>* result, ProbabilityTable* probabilities, bool delete_processed_variants = false);
	/** generates empty UniwueKmers objects for each position (no kmers, only paths). Ownership of vector is transferred to caller. **/
	void compute_empty(std::vector<UniqueKmers*>* result) const;

private:
	KmerCounter* genomic_kmers;
	std::shared_ptr<KmerCounter> read_kmers;
	std::shared_ptr<Graph> variants;
	std::string chromosome;
	size_t kmer_coverage;

	/** compute local coverage in given interval based on unique kmers 
	* @param chromosome chromosome
	* @param var_index variant index
	* @param length how far to go left and right of the variant
	* @returns computed coverage
	**/
	unsigned short compute_local_coverage(size_t var_index, size_t length);
};

#endif // UNIQUEKMERCOMPUTER_HPP
