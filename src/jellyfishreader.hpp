#ifndef JELLYFISHREADER_HPP
#define JELLYFISHREADER_HPP

#include <map>
#include <vector>
#include <string>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
//#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/jellyfish.hpp>
#include <memory>
#include "kmercounter.hpp"

/**
* Reads a kmer counts from the provided Jellyfish .jf file
 **/

/**
struct DataBase {
	private:
		jellyfish::mer_dna_bloom_counter* bloom;
		binary_query* binary;
		bool is_bloom;
		bool is_binary;
	public:
		DataBase() : bloom(nullptr), binary(nullptr), is_bloom(false), is_binary(false) {}
		DataBase(jellyfish::mer_dna_bloom_counter* b) : bloom(b), binary(nullptr), is_bloom(true), is_binary(false) {}
		DataBase(binary_query* b) : bloom(nullptr), binary(b), is_bloom(false), is_binary(true) {}
		size_t check(jellyfish::mer_dna kmer) {
			if (is_bloom) {
				return bloom->check(kmer);
			} else if (is_binary) {
				return binary->check(kmer);
			} else {
				return 0;
			}
		}
		bool next() {
			if (is_bloom) return bloom->next();
			if (is_binary) return binary->next();
			return false;
		}
		size_t val() {
			if (is_bloom) return bloom->val();
			if (is_binary) return bloom->val();
			return 0;
		}
		~DataBase() {
			if (bloom != nullptr) delete bloom;
			if (binary != nullptr) delete binary;
		}
};
**/

class JellyfishReader : public KmerCounter {
public:
	/** 
	* @param readfile name of the FASTQ-files containing reads
	**/
	JellyfishReader(std::string readfile, size_t kmersize);
	
	/** get the abundance of given kmer (string) **/
	size_t getKmerAbundance(std::string kmer);

	/** get the abundance of given kmer (jellyfish kmer) **/
	size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer);

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	size_t computeKmerCoverage(size_t genome_kmers);

	/** computes kmer abundance histogram and returns the three highest peaks **/
	size_t computeHistogram(size_t max_count, bool largest_peak, std::string filename = "");

	~JellyfishReader();

private:
	/** name of input .jf file **/
	std::string filename;
	/** binary map **/
	std::shared_ptr<jellyfish::mapped_file> binary_map;
	/** file header **/
	std::shared_ptr<jellyfish::file_header> header;
	/** jf database **/
	std::shared_ptr<binary_query> db;
	/** infile **/
	std::ifstream ifs;
	
};
#endif // JELLYFISHREADER_HPP
