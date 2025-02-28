#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include <utility>
#include "copynumber.hpp"
#include "kmerpath.hpp"
#include <cereal/types/polymorphic.hpp>
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>


/*
* Represents the set of unique kmers for a variant position.
*/


class UniqueKmers {
public:
	virtual size_t get_variant_position() const = 0;
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	virtual void insert_kmer(unsigned short readcount, std::vector<unsigned short>& allele_ids) = 0;
	/** checks if kmer at index kmer_index is on path path_id **/
	virtual bool kmer_on_path(size_t kmer_index, size_t path_id) const = 0;
	/**  checks if kmer at index kmer_index is on allele allele_id **/
	virtual bool kmer_on_allele(size_t kmer_index, size_t allele_id) const = 0;
	virtual unsigned short get_readcount_of(size_t kmer_index) = 0;
	/** modify kmer count of an already inserted kmer **/
	virtual void update_readcount(size_t kmer_index, unsigned short new_count) = 0;
	/** number of unique kmers **/
	virtual size_t size() const = 0;
	/** return number of paths **/
	virtual unsigned short get_nr_paths() const = 0;
	/** get all paths and alleles covering this position. If only_include, make sure to only output path_ids that are contained in only_include. **/
	virtual void get_path_ids(std::vector<unsigned short>& paths, std::vector<unsigned short>& alleles, std::vector<unsigned short>* only_include = nullptr) = 0;
	/** get all unique alleles covered at this position **/
	virtual void get_allele_ids(std::vector<unsigned short>& a) = 0;
	/** get only those unique alleles which are not undefined **/
	virtual void get_defined_allele_ids(std::vector<unsigned short>& a) = 0;
	/** set the local kmer coverage computed for this position **/
	virtual void set_coverage(unsigned short local_coverage) = 0;
	/** returns the local kmer coverage **/
	virtual unsigned short get_coverage() const = 0;
	/** returns a map which contains the number of unique kmers covering each allele **/
	virtual std::map<unsigned short, int> kmers_on_alleles () const = 0;
	/** returns the number of unique kmers on given allele */
	virtual unsigned short kmers_on_allele(unsigned short allele_id) const = 0;
	/** returns the number of read-supported kmers on given allele **/
	virtual unsigned short present_kmers_on_allele(unsigned short allele_id) const = 0;
	/** returns the fraction of read-supported kmers on given allele **/
	virtual float fraction_present_kmers_on_allele(unsigned short allele_id) const = 0;
	/** check whether allele is undefined **/
	virtual bool is_undefined_allele (unsigned short allele_id) const = 0;
	/** set allele to undefined **/
	virtual void set_undefined_allele (unsigned short allele_id) = 0;
	/** look up allele covered by a path **/
	virtual unsigned short get_allele(unsigned short path_id) const = 0;
	/** update UniqueKmers object by keeping only the paths provided **/
	virtual void update_paths(std::vector<unsigned short>& path_ids) = 0;
	/** print kmer matrix (mainly for debugging) */
	virtual void print_kmer_matrix(std::string chromosome) const = 0;
};


# endif // UNIQUEKMERS_HPP
