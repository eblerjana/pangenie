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
	void add_to_likelihood(unsigned char allele1, unsigned char allele2, long double value);
	/** add allele on haplotype 1 **/
	void add_first_haplotype_allele(unsigned char allele);
	/** add allele on haplotype 2 **/
	void add_second_haplotype_allele(unsigned char allele);
	/** get likelihood of genotype allele1/allele2 (=allele2/allele1) **/
	long double get_genotype_likelihood(unsigned char allele1, unsigned char allele2) const;
	/** get all likelihoods ordered as defined in VCF specification. **/
	std::vector<long double> get_all_likelihoods (size_t nr_alleles) const;
	/** get genotype quality (phred scaled prob that genotype is wrong) **/
	size_t get_genotype_quality (unsigned char allele1, unsigned char allele2) const;
	/** get haplotype **/
	std::pair<unsigned char, unsigned char> get_haplotype() const;
	/** divide all likelihoods by given value **/
	void divide_likelihoods_by(long double value);
	/** get genotype with highest likelihood **/
	std::pair<int, int> get_likeliest_genotype() const;
	/** return number of unique kmers **/
	size_t get_nr_unique_kmers() const;
	/** set number of unique kmers **/
	size_t set_nr_unique_kmers(size_t nr_unique);
	friend std::ostream& operator<<(std::ostream& os, const GenotypingResult& res);
private:
	/** map genotype -> likelihood. genotype alleles are ordered in ascending order **/
	std::map < std::pair<unsigned char,unsigned char>, long double > genotype_to_likelihood;
	unsigned char haplotype_1;
	unsigned char haplotype_2;
	size_t nr_unique_kmers;
};
#endif // GENOTYPINGRESULT_HPP
