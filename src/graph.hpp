#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <memory>
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/vector.hpp>
#include "variant.hpp"
#include "graphbuilder.hpp"
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "uniquekmers.hpp"


template<class T>
std::vector<unsigned char> graph_construct_index(std::vector<T>& alleles, bool reference_added) {
	size_t length = alleles.size();
	unsigned char offset = 0;
	if (reference_added) {
		assert(length > 0);
		length -= 1;
		offset += 1;
	}
	std::vector<unsigned char> index(length);
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), [&](unsigned char a, unsigned char b) { return alleles[a+offset] < alleles[b+offset]; });
	return index;
}

class Graph {
public:

	Graph() = default;
	/**
	* @param fasta_reader FastaReader object containing the reference sequence of the chromosome underlying the graph
	* @param chromosome chromosome the graph represents
	* @param kmer_size kmer size used to build the graph
	* @param reference_added true if reference is used as an additional path in the graph
	**/
	Graph (FastaReader fasta_reader, std::string chromosome, size_t kmer_size, bool reference_added);
	/** returns the kmer size **/
	size_t get_kmer_size() const;
	/** returns the chromosome **/
	std::string get_chromosome() const;
	/** returns the number of variant bubbles in the graph (i.e. after merging clusters) **/
	size_t size() const;
	/** add variant cluster to the graph. Variants are assumed to belong to the same bubble and will thus be merged into one. **/
	void add_variant_cluster(std::vector< std::shared_ptr<Variant> >* cluster, std::vector<std::vector<std::string>>& variant_ids, bool only_defined_ids=false);
	/** return variant bubble at provided index **/
	const Variant& get_variant(size_t index) const;
	/** return FastaReader underlying this object **/
	const FastaReader& get_fasta_reader() const;
	/** write genotyping results into a VCF file **/
	void write_genotypes(std::string filename, const std::vector<GenotypingResult>& genotyping_result, bool write_header, std::string sample, bool ignore_imputed = false);
	/** write phasing results into a VCF file **/
	void write_phasing(std::string filename, const std::vector<GenotypingResult>& genotyping_result, bool write_header, std::string sample, bool ignore_imputed = false);
	/** construct reference sequence left of variant bubble at index **/
	void get_left_overhang(size_t index, size_t length, DnaSequence& result) const;
	/** construct reference sequence right of variant bubble **/
	void get_right_overhang(size_t index, size_t length, DnaSequence& result) const;
	/** deletes a specific variant. This can be used to delete information no longer needed to save space.
	* NOTE: most class functions can no longer be called on an object modified by this function,
	* resulting in an error message.
	*/
	void delete_variant(size_t index);
	/** returns true if at least one variant was deleted by delete_variant function **/
	bool variants_were_deleted() const;

	template<class Archive>
	void save(Archive& archive) const {
		archive(fasta_reader, chromosome, kmer_size, add_reference, variants_deleted, variants, variant_ids);
	}

	template<class Archive>
	void load(Archive& archive) {
		archive(fasta_reader, chromosome, kmer_size, add_reference, variants_deleted, variants, variant_ids);
	}

private:
	/** FastaReader containing reference information **/
	FastaReader fasta_reader;
	/** which chromosome is represented**/
	std::string chromosome;
	/** kmer size **/
	size_t kmer_size;
	/** whether or not reference was added as an additional path **/
	bool add_reference;
	/** indicates whether variants were deleted by delete_variant function **/
	bool variants_deleted;
	/** list of variant bubbles. Individual variants were merged into bubbles. **/
	std::vector<std::shared_ptr<Variant>> variants;
	/** variant IDs in input VCF for all individual variants (i.e. before merging). Needed for output VCF. **/
	std::vector<std::vector<std::string>> variant_ids;

	void insert_ids(std::vector<DnaSequence>& alleles, std::vector<std::string>& variant_ids, bool reference_added);
	std::string get_ids(std::vector<std::string>& alleles, size_t variant_index, bool reference_added);
	friend cereal::access;
};

#endif // GRAPH_HPP
