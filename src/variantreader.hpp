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


template<class T>
std::vector<unsigned char> construct_index(std::vector<T>& alleles, bool reference_added) {
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

class VariantReader {
public:
	VariantReader() = default;
	VariantReader (std::string filename, std::string reference_filename, size_t kmer_size, bool add_reference, std::string sample = "sample");
	~VariantReader();
	/**  writes all path segments (allele sequences + reference sequences in between)
	*    to the given file.
	**/
	size_t get_kmer_size() const;
	void write_path_segments(std::string filename) const;
	void get_chromosomes(std::vector<std::string>* result) const;
	size_t size_of(std::string chromosome) const;
	const Variant& get_variant(std::string chromosome, size_t index) const;
	void open_genotyping_outfile(std::string outfile_name);
	void open_phasing_outfile(std::string outfile_name);
	void write_genotypes_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, bool ignore_imputed = false);
	void write_phasing_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, bool ignore_imputed = false);
	void close_genotyping_outfile();
	void close_phasing_outfile();
	size_t nr_of_paths() const;
	void get_left_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
	void get_right_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
	/** deletes a specific variant. This can be used to delete information no longer needed to save space.
	* NOTE: most class functions can no longer be called on an object modified by this function,
	* resulting in an error message.
	*/
	void delete_variant(std::string chromosome, size_t index);

	template<class Archive>
	void save(Archive& archive) const {
		archive(fasta_reader, kmer_size, nr_paths, nr_variants, add_reference, sample, genotyping_outfile_open, phasing_outfile_open, variants_deleted, variants_per_chromosome, variant_ids);
	}

	template<class Archive>
	void load(Archive& archive) {
		archive(fasta_reader, kmer_size, nr_paths, nr_variants, add_reference, sample, genotyping_outfile_open, phasing_outfile_open, variants_deleted, variants_per_chromosome, variant_ids);
	}

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
	// indicates whether variants were deleted by delete_variant function
	bool variants_deleted;
	std::map< std::string, std::vector< std::shared_ptr<Variant> > > variants_per_chromosome;
	std::map< std::string, std::vector<std::vector<std::string>>> variant_ids;
	void add_variant_cluster(std::string& chromosome, std::vector< std::shared_ptr<Variant> >* cluster);
	void insert_ids(std::string& chromosome, std::vector<DnaSequence>& alleles, std::vector<std::string>& variant_ids, bool reference_added);
	std::string get_ids(std::string chromosome, std::vector<std::string>& alleles, size_t variant_index, bool reference_added);
	friend cereal::access;
};

#endif // VARIANT_READER_HPP
