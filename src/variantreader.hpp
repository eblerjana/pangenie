#ifndef VARIANT_READER_HPP
#define VARIANT_READER_HPP

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
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "uniquekmers.hpp"


class VariantReader {
public:
	VariantReader() = default;
	/** 
	* @param filename name of the input VCF file
	* @param reference_filename name of the reference FASTA file
	* @param result vector to store the constructed Graph segments to (one per chromosome)
	* @param segments_file name of the file to write graph sequences to
	* @param kmer_size kmer size to be used
	* @param add_reference whether or not to add reference as an additional path
	* @param sample name of the sample
	**/
	VariantReader (std::string filename, std::string reference_filename, std::map<std::string, std::shared_ptr<Graph>>& result, std::string segments_file,  size_t kmer_size, bool add_reference, std::string sample = "sample");
	/** return the kmer size **/
	size_t get_kmer_size() const;
	/** get the chromosomes containing variation **/
	void get_chromosomes(std::vector<std::string>* result) const;
	/** total number of paths in the panel **/
	size_t nr_of_paths() const;

private:
	size_t kmer_size;
	size_t nr_variants;
	size_t nr_paths;
	std::vector<string> chromosomes;

	/** reads variants from the VCF and constructs a Graph object for each chromosome. In the Graph, variants closer than kmer_size
	* are merged into a single Variant record that keeps track of the merging (so that in can be undone later) **/
	void construct_graph(std::string filename, std::string reference_filename, std::map<std::string, std::shared_ptr<Graph>>& result, std::string segments_file, std::string sample);
	/** writes a FASTA containing all graph sequences (useful for kmer counting of graph kmers) **/
	void write_path_segments(std::string filename, std::map<std::string, std::shared_ptr<Graph>>& result) const;
	friend cereal::access;
};

#endif // VARIANT_READER_HPP
