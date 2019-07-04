#ifndef GENOTYPINGRESULT_HPP
#define GENOTYPINGRESULT_HPP

#include <map>
#include <utility>
#include <iostream>
#include <vector>

/** Represents the genotyping/phasing result of a position. **/

class GenotypingResult {
public:
	GenotypingResult();
	/** add value to genotype likelihood
	* @param allele1 first genotype allele (arbitrary order)
	* @param allele2 second genotype allele
	* @param value to add to the genotype likelihood
	**/
	void add_to_likelihood(size_t allele1, size_t allele2, long double value);
	/** add allele on haplotype 1 **/
	void add_first_haplotype_allele(size_t allele);
	/** add allele on haplotype 2 **/
	void add_second_haplotype_allele(size_t allele);
	/** get likelihood of genotype allele1/allele2 (=allele2/allele1) **/
	long double get_genotype_likelihood(size_t allele1, size_t allele2) const;
	/** get all likelihoods ordered as defined in VCF specification. **/
	std::vector<long double> get_all_likelihoods (size_t nr_alleles) const;
	/** get genotype quality (phred scaled prob that genotype is wrong) **/
	size_t get_genotype_quality (size_t allele1, size_t allele2) const;
	/** get haplotype **/
	std::pair<size_t, size_t> get_haplotype() const;
	/** divide all likelihoods by given value **/
	void divide_likelihoods_by(long double value);
	/** get genotype with highest likelihood **/
	std::pair<int, int> get_likeliest_genotype() const;
	friend std::ostream& operator<<(std::ostream& os, const GenotypingResult& res);
private:
	/** map genotype -> likelihood. genotype alleles are ordered in ascending order **/
	std::map < std::pair<size_t,size_t>, long double > genotype_to_likelihood;
	size_t haplotype_1;
	size_t haplotype_2;
};
#endif // GENOTYPINGRESULT_HPP
