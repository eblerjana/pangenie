#ifndef GENOTYPINGRESULT_HPP
#define GENOTYPINGRESULT_HPP

#include "uniquekmers.hpp"
#include <map>
#include <utility>
#include <iostream>
#include <vector>
#include <cereal/access.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

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

/** Represents the genotyping/phasing result of a position. **/

class GenotypingResult {
public:
	GenotypingResult();
	/** add value to genotype likelihood
	* @param allele1 first genotype allele (arbitrary order)
	* @param allele2 second genotype allele
	* @param value to add to the genotype likelihood
	**/
	void add_to_likelihood(unsigned short allele1, unsigned short allele2, long double value);
	/** add allele on haplotype 1 **/
	void add_first_haplotype_allele(unsigned short allele);
	/** add allele on haplotype 2 **/
	void add_second_haplotype_allele(unsigned short allele);
	/** get likelihood of genotype allele1/allele2 (=allele2/allele1). **/
	long double get_genotype_likelihood(unsigned short allele1, unsigned short allele2) const;
	/** get all likelihoods ordered as defined in VCF specification. **/
	std::vector<long double> get_all_likelihoods (size_t nr_alleles) const;
	/** get all likelihoods for genotypes containing the given alleles. Likelihoods are normalized so sum up to 1. 
	NOTE: haplotype alleles are set only if they occur in the list of given alleles. Otherwise (i.e. if undefined), they are 0.**/
	GenotypingResult get_specific_likelihoods (std::vector<unsigned short>& alleles) const;
	/** get genotype quality (phred scaled prob that genotype is wrong) **/
	size_t get_genotype_quality (unsigned short allele1, unsigned short allele2) const;
	/** get haplotype **/
	std::pair<unsigned short, unsigned short> get_haplotype() const;
	/** divide all likelihoods by given value **/
	void divide_likelihoods_by(long double value);
	/** get genotype with highest likelihood **/
	std::pair<int, int> get_likeliest_genotype() const;
	friend std::ostream& operator<<(std::ostream& os, const GenotypingResult& res);
	/** add the given likelihoods. Haplotypes, kmer_counts and coverage are not changed,
	 ** only the likelihoods are modified and not normalized.
	 **/
	void combine(GenotypingResult& likelihoods);
	void normalize();
	void set_unique_kmers(unsigned short nr_unique_kmers);
	void set_coverage(unsigned short coverage);
	unsigned short nr_unique_kmers() const;
	unsigned short coverage() const;
	/** returns true if there are no likelihoods stored **/
	bool contains_no_likelihoods() const;

	/** provide access to stored likelihoods **/
	const std::map < std::pair<unsigned short,unsigned short>, long double >& get_stored_likelihoods() const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(genotype_to_likelihood, haplotype_1, haplotype_2, local_coverage, unique_kmers);
	}

private:
	/** map genotype -> likelihood. genotype alleles are ordered in ascending order **/
	std::map < std::pair<unsigned short,unsigned short>, long double > genotype_to_likelihood;
	unsigned short haplotype_1;
	unsigned short haplotype_2;
	unsigned short local_coverage;
	unsigned short unique_kmers;
	friend cereal::access;
};
#endif // GENOTYPINGRESULT_HPP
