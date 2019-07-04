#ifndef FASTAREADER_HPP
#define FASTAREADER_HPP

#include <string>
#include <map>
#include "dnasequence.hpp"

/**
* Represents FASTA-sequence.
**/

class FastaReader {
public:
	/** 
	* @param filename name of the FASTA-file.
	**/
	FastaReader(std::string filename);
	~FastaReader();
	/** check if sequence with given name exists in file. **/
	bool contains_name(std::string name) const;
	/** get length of sequence with name **/
	size_t get_size_of(std::string name) const;
	/** compute total numbers of kmers in the sequences **/
	size_t get_total_kmers(size_t kmer_size) const;
	void get_subsequence(std::string name, size_t start, size_t end, std::string& result) const;
private:
	void parse_file(std::string filename);
	std::map<std::string, DnaSequence*> name_to_sequence;
};

#endif // FASTAREADER_HPP
