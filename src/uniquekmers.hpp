#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <vector>
#include <string>
#include <map>
#include <utility>
#include "copynumber.hpp"
#include "kmerpath.hpp"
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>


/*
* Represents the set of unique kmers for a variant position.
*/

// serialization of std::pair, code taken from: https://github.com/USCiLab/cereal/issues/547
namespace cereal
{
    template<class Archive, class F, class S>
    void save(Archive& ar, const std::pair<F, S>& pair)
    {
        ar(pair.first, pair.second);
    }

    template<class Archive, class F, class S>
    void load(Archive& ar, std::pair<F, S>& pair)
    {
        ar(pair.first, pair.second);
    }

    template <class Archive, class F, class S> 
    struct specialize<Archive, std::pair<F, S>, cereal::specialization::non_member_load_save> {};
}


class UniqueKmers {
public:
	/**
	* @param variant_position genomic variant position
	* @param alleles defines which path (= index) covers each allele (= alleles[index])
	**/
	UniqueKmers() = default;
	UniqueKmers(size_t variant_position, std::vector<unsigned char>& alleles);
	/** returns the variant position **/
	size_t get_variant_position();
	/** insert a kmer
	* @param cn copy number probabilities of kmer
	* @param allele_ids on which alleles this kmer occurs
	**/
	void insert_kmer(unsigned short readcount, std::vector<unsigned char>& allele_ids);
	/** checks if kmer at index kmer_index is on path path_id **/
	bool kmer_on_path(size_t kmer_index, size_t path_id) const;
	unsigned short get_readcount_of(size_t kmer_index);
	/** modify kmer count of an already inserted kmer **/
	void update_readcount(size_t kmer_index, unsigned short new_count);
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
	/** look up allele covered by a path **/
	unsigned char get_allele(unsigned short path_id) const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(variant_pos, current_index, kmer_to_count, alleles, path_to_allele, local_coverage);
	}

private:
	size_t variant_pos;
	size_t current_index;
	std::vector<unsigned short> kmer_to_count;
	// stores kmers of each allele and whether the allele is undefined
	std::map<unsigned char, std::pair<KmerPath, bool>> alleles;
	// defines which alleles are carried by each path (=index)
	std::vector<unsigned char> path_to_allele;
	unsigned short local_coverage;
	friend class EmissionProbabilityComputer;
	friend cereal::access;
};
# endif // UNIQUEKMERS_HPP
