#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
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
	UniqueKmers(size_t variant_id, size_t variant_position);
	size_t get_variant_index();
	size_t get_variant_position();
	/** insert empty allele (no kmers) **/
	void insert_empty_allele(unsigned char allele_id, bool is_undefined = false);
	/** insert a path covering the given allele **/
	void insert_path(size_t path_id, unsigned char allele_id);
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(CopyNumber cn, std::vector<unsigned char>& allele_ids);
	/** checks if kmer at index kmer_index is on path path_id **/
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
	/** get only those unique alleles which are not undefined **/
	void get_defined_allele_ids(std::vector<unsigned char>& a);
	friend std::ostream& operator<< (std::ostream& stream, const UniqueKmers& uk);
	/** set the local kmer coverage computed for this position **/
	void set_coverage(double local_coverage);
	/** returns the local kmer coverage **/
	double get_coverage() const;
	/** returns a map which contains the number of unique kmers covering each allele **/
	std::map<unsigned char, int> kmers_on_alleles () const;
	/** check whether allele is undefined **/
	bool is_undefined_allele (unsigned char allele_id) const;
	/** set allele to undefined **/
	void set_undefined_allele (unsigned char allele_id);

private:
	size_t variant_id;
	size_t variant_pos;
	size_t current_index;
	std::vector<CopyNumber> kmer_to_copynumber;
	std::map<unsigned char, KmerPath> alleles;
	/** keep track of whether an allele is undefined **/
	std::map<unsigned char, bool> is_undefined;
	std::map<size_t, unsigned char> path_to_allele;
	CopyNumberAssignment combine_paths(unsigned char allele_id1, unsigned char allele_id2);
	double local_coverage;
	friend class EmissionProbabilityComputer;
	
};
# endif // UNIQUEKMERS_HPP
