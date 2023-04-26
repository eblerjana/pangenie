#ifndef GRAPH_BUILDER_HPP
#define GRAPH_BUILDER_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <memory>
#include "graph.hpp"
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "uniquekmers.hpp"


class GraphBuilder {
public:
	GraphBuilder() = default;
	/** 
	* @param filename name of the input VCF file
	* @param reference_filename name of the reference FASTA file
	* @param result vector to store the constructed Graph segments to (one per chromosome)
	* @param segments_file name of the file to write graph sequences to
	* @param kmer_size kmer size to be used
	* @param add_reference whether or not to add reference as an additional path
	**/
	GraphBuilder (std::string filename, std::string reference_filename, std::map<std::string, std::shared_ptr<Graph>>& result, std::string segments_file,  size_t kmer_size, bool add_reference);
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
	std::vector<std::string> chromosomes;

	/** reads variants from the VCF and constructs a Graph object for each chromosome. In the Graph, variants closer than kmer_size
	* are merged into a single Variant record that keeps track of the merging (so that in can be undone later) **/
	void construct_graph(std::string filename, FastaReader* fasta_reader, std::map<std::string, std::shared_ptr<Graph>>& result, bool add_reference);
	/** writes a FASTA containing all graph sequences (useful for kmer counting of graph kmers) **/
	void write_path_segments(std::string filename, FastaReader* fasta_reader, std::map<std::string, std::shared_ptr<Graph>>& result) const;
};

#endif // GRAPH_BUILDER_HPP
