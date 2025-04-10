#ifndef STEPWISEUNIQUEKMERCOMPUTER_HPP
#define STEPWISEUNIQUEKMERCOMPUTER_HPP

#include <vector>
#include <string>
#include <memory>
#include <zlib.h>
#include "kmercounter.hpp"
#include "graph.hpp"
#include "multiallelicuniquekmers.hpp"
#include "biallelicuniquekmers.hpp"
#include "probabilitytable.hpp"

class StepwiseUniqueKmerComputer {
public:
	/** 
	* @param genomic_kmers genomic kmer counts
	* @param variants
	**/
	StepwiseUniqueKmerComputer (KmerCounter* genomic_kmers, std::shared_ptr<Graph> variants);
	/** generates UniqueKmers object for each position, ownership of vector is transferred to the caller.
	* @param result	UniqueKmer objects will be stored here
	* @param filename name of file to write kmer information to
	* @param delete_processed_variants if set to true, the Graph will be motified. Processed Variant objects will be deleted. Use with care!
	 **/
	void compute_unique_kmers(std::vector<std::shared_ptr<UniqueKmers>>* result, std::string filename,  bool delete_processed_variants = false);
	/** generates UniqueKmers object for each position, ownership of vector is transferred to the caller. Writes output data in FASTA format.
	* @param result	UniqueKmer objects will be stored here
	* @param filename name of file to write kmer information to
	* @param delete_processed_variants if set to true, the Graph will be motified. Processed Variant objects will be deleted. Use with care!
	 **/

	/** generates empty UniwueKmers objects for each position (no kmers, only paths). Ownership of vector is transferred to caller. **/
	void compute_empty(std::vector<std::shared_ptr<UniqueKmers>>* result) const;

private:
	KmerCounter* genomic_kmers;
	std::shared_ptr<Graph> variants;
	std::string chromosome;
	void determine_unique_flanking_kmers(size_t var_index, size_t length, std::vector<std::string>& result);

	/**
	* @param occurances maps unique kmers to a list of alleles they occur on
	* @returns map of alleles to a list of selected unique kmers
	**/
	std::map<unsigned short, std::vector<jellyfish::mer_dna>> select_kmers(const Variant* variant, std::map <jellyfish::mer_dna, std::vector<unsigned short>>& occurences, bool is_biallelic);
};

#endif // STEPWISEUNIQUEKMERCOMPUTER_HPP
