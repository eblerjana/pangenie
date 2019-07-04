#ifndef DNASEQUENCE_HPP
#define DNASEQUENCE_HPP

#include <string>
#include <vector>

/** 
* Represents a DNA sequence.
**/

class DnaSequence {
public:
	DnaSequence();
	DnaSequence(std::string& sequence);
	/** append sequence to the end of DNA. **/
	void append(std::string& sequence);
	/** get base at index position. **/
	char operator[](size_t position) const;
	size_t length() const;
	/** get subtring 
	* @param start, end start and end of the subsequence
	* @param result resulting subtring
	**/
	void substr(size_t start, size_t end, std::string& result) const;

private:
	/** store 2 bases per char (using 4 bits for each) **/
	std::vector<unsigned char> sequence;
	/** length of the sequence **/
	size_t size;
};

#endif // DNASEQUENCE_HPP
