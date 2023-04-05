#ifndef DNASEQUENCE_HPP
#define DNASEQUENCE_HPP

#include <string>
#include <vector>
#include <memory>
#include <cereal/access.hpp>
#include <cereal/types/vector.hpp>

/** 
* Represents a DNA sequence.
**/

class DnaSequence {
public:
	DnaSequence();
	DnaSequence(std::string& sequence);
	/** append string sequence to the end of DNA. **/
	void append(std::string& sequence);
	/** append DnaSequence to the end of DNA **/
	void append(DnaSequence sequence);
	/** reverse sequence **/
	void reverse();
	/** compute reverse complement **/
	void reverse_complement();
	/** get base at index position. **/
	char operator[](size_t position) const;
	DnaSequence base_at(size_t position) const;
	size_t size() const;
	/** get subsequence 
	* @param start, end start and end of the subsequence
	* @param result resulting DnaSequence
	**/
	void substr(size_t start, size_t end, DnaSequence& result) const;
	/** get subsequence a string **/
	void substr(size_t start, size_t end, std::string& result) const;
	/** convert DnaSequence to string **/
	std::string to_string() const;
	/** clear sequence **/
	void clear();
	bool operator<(const DnaSequence& dna) const;
	/** comparision operators **/
	friend bool operator==(const DnaSequence& dna1, const DnaSequence& dna2);
	friend bool operator!=(const DnaSequence& dna1, const DnaSequence& dna2);
	/** returns true if the sequence contains at least one undefined base (different from A,T,C,G,a,t,c,g) **/
	bool contains_undefined() const;

	template<class Archive>
	void serialize(Archive& archive) {
		archive(sequence, even_length, is_undefined);
	}

private:
	/** store 2 bases per char (using 4 bits for each) **/
	std::vector<unsigned char> sequence;
	/** true if the length is even **/
	bool even_length;
	/** is true if sequence contains any undefined bases (N,n) **/
	bool is_undefined;
	friend cereal::access;
};

#endif // DNASEQUENCE_HPP
