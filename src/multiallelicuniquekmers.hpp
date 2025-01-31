#ifndef MULTIALLELICUNIQUEKMERS_HPP
#define MULTIALLELICUNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include <utility>
#include "uniquekmers.hpp"
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


class MultiallelicUniqueKmers : public UniqueKmers {
public:
	/**
	* @param variant_position genomic variant position
	* @param alleles defines which path (= index) covers each allele (= alleles[index])
	**/
	MultiallelicUniqueKmers() = default;
	MultiallelicUniqueKmers(size_t variant_position, std::vector<unsigned short>& alleles);

	size_t get_variant_position() const;
	void set_coverage(unsigned short local_coverage); 
	unsigned short get_coverage() const;
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(unsigned short readcount, std::vector<unsigned short>& allele_ids);
	/** checks if kmer at index kmer_index is on path path_id **/
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	/** checks if kmer at index kmer_index is on allele allele_id 
	* Note: optimized for speed, no boundary checks performed.
	**/
	bool kmer_on_allele(size_t kmer_index, size_t allele_id) const;
	unsigned short get_readcount_of(size_t kmer_index);
	/** modify kmer count of an already inserted kmer **/
	void update_readcount(size_t kmer_index, unsigned short new_count);
	/** number of unique kmers **/
	size_t size() const;
	/** return number of paths **/
	unsigned short get_nr_paths() const;
	/** get all paths and alleles covering this position. If only_include, make sure to only output path_ids that are contained in only_include. **/
	void get_path_ids(std::vector<unsigned short>& paths, std::vector<unsigned short>& alleles, std::vector<unsigned short>* only_include = nullptr);
	/** get all unique alleles covered at this position **/
	void get_allele_ids(std::vector<unsigned short>& a);
	/** get only those unique alleles which are not undefined **/
	void get_defined_allele_ids(std::vector<unsigned short>& a);
	friend std::ostream& operator<< (std::ostream& stream, const MultiallelicUniqueKmers& uk);
	/** returns a map which contains the number of unique kmers covering each allele **/
	std::map<unsigned short, int> kmers_on_alleles () const;
	/** returns the number of unique kmers on given allele */
	unsigned short kmers_on_allele(unsigned short allele_id) const;
	/** returns the number of read-supported kmers on given allele **/
	unsigned short present_kmers_on_allele(unsigned short allele_id) const;
	/** returns the fraction of read-supported kmers on given allele **/
	float fraction_present_kmers_on_allele(unsigned short allele_id) const;
	/** check whether allele is undefined **/
	bool is_undefined_allele (unsigned short allele_id) const;
	/** set allele to undefined **/
	void set_undefined_allele (unsigned short allele_id);
	/** look up allele covered by a path **/
	unsigned short get_allele(unsigned short path_id) const;
	/** update MultiallelicUniqueKmers object by keeping only the paths provided **/
	void update_paths(std::vector<unsigned short>& path_ids);
	/** print kmer matrix (mainly for debugging) */
	void print_kmer_matrix(std::string chromosome) const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(variant_pos, local_coverage, current_index, kmer_to_count, alleles, path_to_allele);
	}

private:
	size_t variant_pos;
	float local_coverage;
	size_t current_index;
	std::vector<unsigned short> kmer_to_count;
	// stores kmers of each allele and whether the allele is undefined
	std::map<unsigned short, AlleleInfo> alleles;
	// defines which alleles are carried by each path (=index)
	std::vector<unsigned short> path_to_allele;
	friend class HaplotypeSampler;
	friend cereal::access;
};

#include <cereal/archives/binary.hpp>
#include <cereal/details/static_object.hpp>
CEREAL_REGISTER_TYPE(MultiallelicUniqueKmers);
CEREAL_REGISTER_POLYMORPHIC_RELATION(UniqueKmers, MultiallelicUniqueKmers)


# endif // MULTIALLELICUNIQUEKMERS_HPP
