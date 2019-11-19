#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include "copynumber.hpp"
#include "kmerpath.hpp"
#include "dynamicbitset.hpp"

/*
* Represents the set of unique kmers for a variant position.
*/


class UniqueKmers {
public:
	/**
	* @param variant_id variant identifier
	* @param variant_position genomic variant position
	**/
	UniqueKmers(size_t variant_id, size_t variant_position);
	size_t get_variant_index();
	size_t get_variant_position();
	/** insert empty allele (no kmers) **/
	void insert_empty_allele(unsigned char allele_id);
	/** insert a path covering the given allele **/
	void insert_path(size_t path_id, unsigned char allele_id);
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(CopyNumber cn, std::vector<unsigned char>& allele_ids);
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	CopyNumber get_copynumber_of(size_t kmer_index);
	/** number of unique kmers **/
	size_t size() const;
	/** return number of paths **/
	size_t get_nr_paths() const;
	/** get all paths and alleles covering this position **/
	void get_path_ids(std::vector<size_t>& paths, std::vector<unsigned char>& alleles);
	/** get all unique alleles covered at this position **/
	void get_allele_ids(std::vector<unsigned char>& a);
	/** get a bitset that encodes all paths that cover the allele allele_id **/
	DynamicBitset get_paths_of_allele(unsigned char allele_id) const;
	friend std::ostream& operator<< (std::ostream& stream, const UniqueKmers& uk);

private:
	size_t variant_id;
	size_t variant_pos;
	size_t current_index;
	size_t nr_paths;
	size_t max_path;
	std::vector<CopyNumber> kmer_to_copynumber;
	std::map<unsigned char, KmerPath> alleles;
//	std::map<size_t, unsigned char> path_to_allele;
	std::map<unsigned char, DynamicBitset> allele_to_paths;

	CopyNumberAssignment combine_paths(unsigned char allele_id1, unsigned char allele_id2);
	friend class EmissionProbabilityComputer;
};
# endif // UNIQUEKMERS_HPP
