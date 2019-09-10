#ifndef JELLYFISHCOUNTER_HPP
#define JELLYFISHCOUNTER_HPP

#include <map>
#include <vector>
#include <string>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>
#include "kmercounter.hpp"

/**
* Counts Kmers in DNA-sequences (given in FASTQ-format) using jellyfish.
 **/


// using example from: Jellyfish-2/examples/jf_count_dump/jf_count_dump.cc
typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>	mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>>	sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna>	mer_iterator_type;

class mer_counter : public jellyfish::thread_exec {
	mer_hash_type& mer_hash_;
	jellyfish::stream_manager<char**> streams_;
	sequence_parser_type parser_;
  	const bool canonical_;

public:
	mer_counter(int nb_threads, mer_hash_type& mer_hash,
	char** file_begin, char** file_end,
	bool canonical)
	: mer_hash_(mer_hash)
	, streams_(file_begin, file_end)
	, parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
	, canonical_(canonical)
{ }

	virtual void start(int thid) {
	mer_iterator_type mers(parser_, canonical_);

	for( ; mers; ++mers)
		mer_hash_.add(*mers, 1);
	mer_hash_.done();
	}
};


class JellyfishCounter : public KmerCounter {
public:
	/** 
	* @param readfile name of the FASTQ-files containing reads
	* @param *params parameters for GATB-Kmercounter
	* @param name of the output file
	**/
	JellyfishCounter(std::string readfile, size_t kmer_size, size_t nr_threads = 1);

	~JellyfishCounter();
	
	/** get the abundance of given kmer (string) **/
	size_t getKmerAbundance(std::string kmer);

	/** get the abundance of given kmer (jellyfish kmer) **/
	size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer);

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	size_t computeKmerCoverage(size_t genome_kmers);

	/** computes kmer abundance histogram and returns the three highest peaks **/
	size_t computeHistogram(size_t max_count, std::string filename = "");

private:
	mer_hash_type* jellyfish_hash;
};
#endif // JELLYFISHCOUNTER_HPP
