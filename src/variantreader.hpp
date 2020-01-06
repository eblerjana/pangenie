#ifndef VARIANT_READER_HPP
#define VARIANT_READER_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"

class VariantReader {
public:
	VariantReader (std::string filename, std::string reference_filename, size_t kmer_size, bool add_reference, std::string sample = "sample");
	/**  writes all path segments (allele sequences + reference sequences in between)
	*    to the given file.
	**/
	size_t get_kmer_size() const;
	void write_path_segments(std::string filename) const;
	void get_chromosomes(std::vector<std::string>* result) const;
	size_t size_of(std::string chromosome) const;
	const Variant& get_variant(std::string chromosome, size_t index) const;
	const std::vector<Variant>& get_variants_on_chromosome(std::string chromosome) const;
	void open_genotyping_outfile(std::string outfile_name);
	void open_phasing_outfile(std::string outfile_name);
	void write_genotypes_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, bool ignore_imputed = false);
	void write_phasing_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, bool ignore_imputed = false);
	void close_genotyping_outfile();
	void close_phasing_outfile();
	size_t nr_of_genomic_kmers() const;
	void get_left_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
	void get_right_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;

private:
	FastaReader fasta_reader;
	size_t kmer_size;
	size_t nr_paths;
	size_t nr_variants;
	bool add_reference;
	std::string sample;
	std::ofstream genotyping_outfile;
	std::ofstream phasing_outfile;
	bool genotyping_outfile_open;
	bool phasing_outfile_open;
	std::map< std::string, std::vector<Variant> > variants_per_chromosome;
	void add_variant_cluster(std::string& chromosome, std::vector<Variant>* cluster);
};

#endif // VARIANT_READER_HPP
