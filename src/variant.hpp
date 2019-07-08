#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <vector>
#include <string>
#include <iostream>
#include "genotypingresult.hpp"

/** 
* Represents a variant.
**/

class Variant {
public:
	/** 
	* @param left_flank, right_flank sequences left and right of the variant
	* @param start_position, end_position coordinates
	* @param alleles list of alleles (first one is used as reference allele)
	* @param paths vector containing the allele each path covers (i-th path covers allele at paths[i])
	**/
	Variant(std::string left_flank, std::string right_flank, std::string chromosome, size_t start_position, size_t end_position, std::vector<std::string> alleles, std::vector<size_t> paths);
	void add_flanking_sequence();
	void remove_flanking_sequence();
	void combine_variants (Variant const &v2);
	void separate_variants (std::vector<Variant>* resulting_variants, const GenotypingResult* input_genotyping = nullptr, std::vector<GenotypingResult>* resulting_genotyping = nullptr) const;
	size_t nr_of_alleles() const;
	size_t nr_of_paths() const;
	std::string get_allele_sequence(size_t index) const;
	size_t get_start_position() const;
	size_t get_end_position() const;
	std::string get_chromosome() const;
	bool allele_on_path(size_t allele_index, size_t path_index) const;
	size_t get_allele_on_path(size_t path_index) const;
	void get_paths_of_allele(size_t allele_index, std::vector<size_t>& result) const;
	bool is_combined() const;
	friend std::ostream& operator<<(std::ostream& os, const Variant& var);
	friend bool operator==(const Variant& v1, const Variant& v2);
	friend bool operator!=(const Variant& v1, const Variant& v2);
private:
	std::string left_flank;
	std::string right_flank;
	std::string chromosome;
	std::vector<size_t> start_positions;
	std::vector<size_t> end_positions;
	std::vector<std::vector<std::string>> alleles;
	/** alleles not covered by any paths **/
	std::vector<std::vector<std::string>> uncovered_alleles;
	std::vector<size_t> paths;
	bool flanks_added;
	
};

#endif //VARIANT_HPP
