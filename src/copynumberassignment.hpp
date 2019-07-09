#ifndef COPYNUMBERASSIGNMENT_HPP
#define COPYNUMBERASSIGNMENT_HPP

#include <vector>
#include <iostream>

/**
* Represents an assignment of copy numbers to a sequence of kmers.
**/

class CopyNumberAssignment {
public:
	CopyNumberAssignment();
	CopyNumberAssignment(std::vector<unsigned int> kmers);
	/** assign kmer at index copy number value. 
	* @param index kmer index
	* @param value copy number (only values 0,1,2 are valid)
	**/
	void set_position(size_t index, unsigned int value);
	/** return assigned copy number at index. **/
	unsigned int get_position(size_t index) const;
	friend std::ostream& operator<< (std::ostream& stream, const CopyNumberAssignment& cna);
	std::string convert_to_string() const;

protected:
	/** use one unsigned int to store the assignments of 20 kmers (base 2). **/
	std::vector<unsigned int> kmers;
};

#endif // COPYNUMBERASSIGNMENT_HPP
