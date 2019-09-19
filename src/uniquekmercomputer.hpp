#ifndef UNIQUEKMERCOMPUTER_HPP
#define UNIQUEKMERCOMPUTER_HPP

#include <vector>
#include <string>
#include "kmercounter.hpp"
#include "variantreader.hpp"
#include "uniquekmers.hpp"
#include "probabilitycomputer.hpp"

class UniqueKmerComputer {
public:
	/** 
	* @param genomic_kmers genomic kmer counts
	* @param read_kmers read kmer counts
	* @param variants 
	* @param kmer_coverage needed to compute kmer copy number probabilities
	**/
	UniqueKmerComputer (KmerCounter* genomic_kmers, KmerCounter* read_kmers, VariantReader* variants, std::string chromosome, size_t kmer_coverage);
	/** generates UniqueKmers object for each position, ownership of vector is transferred to the caller. **/
	void compute_unique_kmers(std::vector<UniqueKmers*>* result, long double regularization_const = 0.0L) const;
	/** generates empty UniwueKmers objects for each position (no kmers, only paths). Ownership of vector is transferred to caller. **/
	void compute_empty(std::vector<UniqueKmers*>* result) const;

private:
	KmerCounter* genomic_kmers;
	KmerCounter* read_kmers;
	VariantReader* variants;
	std::string chromosome;
	ProbabilityComputer probability;
	size_t kmer_coverage;
};

#endif // UNIQUEKMERCOMPUTER_HPP
