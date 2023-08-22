#ifndef VARIANT_HPP
#define VARIANT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <cereal/access.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include "genotypingresult.hpp"
#include "dnasequence.hpp"
#include "uniquekmers.hpp"

/** 
* Represents a variant.
**/

struct VariantStats {
	size_t nr_unique_kmers;
	std::map<unsigned char, int> kmer_counts;
	unsigned short coverage;
};

class Graph;

class Variant {
public:
	/** 
	* @param left_flank, right_flank sequences left and right of the variant
	* @param chromosome which chromosome the variant is located on
	* @param start_position, end_position coordinates
	* @param alleles list of alleles (first one is reference allele)
	* @param paths vector containing the allele each path covers (i-th path covers allele at paths[i])
	* @param variant_id ID of the variant (ID column of the VCF)
	*
	* Currently, the largest number of alleles and paths supported is 256. This is because unsigned chars are used to store allele ids. For merged variants, the number of alleles can get as
	* high as the number of paths (if every path carries a different allele), therefore the limit on the number of paths is also 256 to avoid overflows.
	**/
	Variant() = default;
	Variant(std::string left_flank, std::string right_flank, std::string chromosome, size_t start_position, size_t end_position, std::vector<std::string> alleles, std::vector<unsigned char> paths); //, std::string variant_id = ".");
	Variant(DnaSequence& left_flank, DnaSequence& right_flank, std::string chromosome, size_t start_position, size_t end_position, std::vector<DnaSequence>& alleles, std::vector<unsigned char>& paths); //, std::string variant_id = ".");
	/** add flanking sequences left and right of variant **/
	void add_flanking_sequence();
	/** remove flanking sequences left and right of variant **/
	void remove_flanking_sequence();
	/** combine variants into a multi-allelic variant **/
	void combine_variants (Variant const &v2);
	/** separate variants that have been combined **/
	void separate_variants (std::vector<Variant>* resulting_variants, const GenotypingResult* input_genotyping = nullptr, std::vector<GenotypingResult>* resulting_genotyping = nullptr) const;
	/** total number of alleles of the variant **/
	size_t nr_of_alleles() const;
	/** total number of paths covering the variant **/
	size_t nr_of_paths() const;
	/** return allele sequence as string **/
	std::string get_allele_string(size_t index) const;
	/** return allele sequence as DnaSequence **/
	DnaSequence get_allele_sequence(size_t index) const;
	/** get start position of the variant **/
	size_t get_start_position() const;
	/** get end position of the variant **/
	size_t get_end_position() const;
	/** get chromosome **/
	std::string get_chromosome() const;
	/** check if given allele is covered by the given path **/
	bool allele_on_path(unsigned char allele_index, size_t path_index) const;
	/** get index of allele located on given path **/
	unsigned char get_allele_on_path(size_t path_index) const;
	/** get paths that cover a given allele **/
	void get_paths_of_allele(unsigned char allele_index, std::vector<size_t>& result) const;
	/** check if this is a combined variant **/
	bool is_combined() const;
	friend std::ostream& operator<<(std::ostream& os, const Variant& var);
	friend bool operator==(const Variant& v1, const Variant& v2);
	friend bool operator!=(const Variant& v1, const Variant& v2);
	/** compute allele frequency of the given allele **/
	float allele_frequency(unsigned char allele_index, bool ignore_ref_path = false) const;
	/** return variant ID **/
	std::string get_id() const;
	/** check whether the given allele is undefined **/
	bool is_undefined_allele(size_t allele_id) const;
	/** return number of paths with missing alleles **/
	size_t nr_missing_alleles() const;
	/** determine statistics for each individual variant **/
	void variant_statistics (std::shared_ptr<UniqueKmers> unique_kmers, std::vector<VariantStats>& result) const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(left_flank, right_flank, inner_flanks, chromosome, start_position, allele_sequences, allele_combinations, uncovered_alleles, paths, flanks_added);
	}

private:
	// flanking sequence at left end
	DnaSequence left_flank;
	// flanking sequence at right end
	DnaSequence right_flank;
	// sequences between two combined alleles
	std::vector<DnaSequence> inner_flanks;
	// chromosome
	std::string chromosome;
	// starting position of variant
	size_t start_position;
	// IDs of individual variants
//	std::vector<std::string> variant_ids;
	// allele sequences of individual variants
	std::vector<std::vector<DnaSequence>> allele_sequences;
	// combined alleles
	std::vector<std::vector<unsigned char>> allele_combinations;
	// alleles not covered by any path
	std::vector<std::vector<unsigned char>> uncovered_alleles;
	std::vector<unsigned char> paths;
	bool flanks_added;
	void set_values(size_t end_position);
	friend cereal::access;
	friend Graph;
};

#endif //VARIANT_HPP
