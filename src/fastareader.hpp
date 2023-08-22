#ifndef FASTAREADER_HPP
#define FASTAREADER_HPP

#include <string>
#include <map>
#include <memory>
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/memory.hpp>
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
	FastaReader() = default;
	/** check if sequence with given name exists in file. **/
	bool contains_name(std::string name) const;
	/** get length of sequence with name **/
	size_t get_size_of(std::string name) const;
	/** get names of seqences present in the file **/
	void get_sequence_names(std::vector<std::string>& names) const;
	/** compute total numbers of kmers in the sequences **/
	size_t get_total_kmers(size_t kmer_size) const;
	/** get a subsequence **/
	void get_subsequence(std::string name, size_t start, size_t end, std::string& result) const;
	void get_subsequence(std::string name, size_t start, size_t end, DnaSequence& result) const;
	/** creates a new FastaReader object containing the sequence of the given sequence name. 
	* This sequence will no longer be stored in the current object (i.e. it will be moved to the
	* resulting FastaReader)
	**/
	FastaReader extract_name(std::string chromosome);

	template<class Archive>
	void serialize(Archive& archive) {
		archive(name_to_sequence);
	}

private:
	void parse_file(std::string filename);
	std::map<std::string, std::shared_ptr<DnaSequence>> name_to_sequence;
	friend cereal::access;
};

#endif // FASTAREADER_HPP
