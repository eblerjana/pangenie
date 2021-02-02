#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include <utility>
#include "copynumber.hpp"
#include "kmerpath.hpp"

/*
* Represents the set of unique kmers for a variant position.
*/


class UniqueKmers {
public:
	/**
	* @param variant_id variant identifier
	* @param variant_position genomic variant position
	**/
	UniqueKmers(size_t variant_position);
	size_t get_variant_position();
	/** insert empty allele (no kmers) **/
	void insert_empty_allele(unsigned char allele_id, bool is_undefined = false);
	/** insert a path covering the given allele **/
	void insert_path(unsigned short path_id, unsigned char allele_id);
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(unsigned short readcount, std::vector<unsigned char>& allele_ids);
	/** checks if kmer at index kmer_index is on path path_id **/
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	unsigned short get_readcount_of(size_t kmer_index);
	/** number of unique kmers **/
	size_t size() const;
	/** return number of paths **/
	unsigned short get_nr_paths() const;
	/** get all paths and alleles covering this position. If only_include, make sure to only output path_ids that are contained in only_include. **/
	void get_path_ids(std::vector<unsigned short>& paths, std::vector<unsigned char>& alleles, std::vector<unsigned short>* only_include = nullptr);
	/** get all unique alleles covered at this position **/
	void get_allele_ids(std::vector<unsigned char>& a);
	/** get only those unique alleles which are not undefined **/
	void get_defined_allele_ids(std::vector<unsigned char>& a);
	friend std::ostream& operator<< (std::ostream& stream, const UniqueKmers& uk);
	/** set the local kmer coverage computed for this position **/
	void set_coverage(unsigned short local_coverage);
	/** returns the local kmer coverage **/
	unsigned short get_coverage() const;
	/** returns a map which contains the number of unique kmers covering each allele **/
	std::map<unsigned char, int> kmers_on_alleles () const;
	/** check whether allele is undefined **/
	bool is_undefined_allele (unsigned char allele_id) const;
	/** set allele to undefined **/
	void set_undefined_allele (unsigned char allele_id);

private:
	size_t variant_pos;
	size_t current_index;
	std::vector<unsigned short> kmer_to_count;
	// stores kmers of each allele and whether the allele is undefined
	std::map<unsigned char, std::pair<KmerPath, bool>> alleles;
	std::map<unsigned short, unsigned char> path_to_allele;
	unsigned short local_coverage;
	friend class EmissionProbabilityComputer;
	
};
# endif // UNIQUEKMERS_HPP
